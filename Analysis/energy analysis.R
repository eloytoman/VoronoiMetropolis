library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)


energy_analysis <- function(points){
  
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Bending_energy = double(Lay),
                         Total_energy=double(Lay))
  
  for (i in 1:Lay) {
    tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
    tesener[i,c(2,3,4,6)] <- c(elener,
                               tenener,
                               contener,
                               sum((areas-1)^2+(gam/2)*(perims^2)+
                                     lamad*perims)/Lay)
  
  }
  
  tesener$Bending_energy[c(1,Lay)]<-0
  
  for (i in 2:(Lay-1)) {
    
    angles <- numeric(n)
    
    for (j in 1:n) {
      ptcentral <- c(points[[i]]$y[j],
                     rad[[i]]*cos((1/rad[[i]])*points[[i]]$x[j]),
                     rad[[i]]*sin((1/rad[[i]])*points[[i]]$x[j]))
      
      ptinf <- c(points[[i-1]]$y[j],
                 rad[[i-1]]*cos((1/rad[[i-1]])*points[[i-1]]$x[j]),
                 rad[[i-1]]*sin((1/rad[[i-1]])*points[[i-1]]$x[j]))
      
      ptsup <- c(points[[i+1]]$y[j],
                 rad[[i+1]]*cos((1/rad[[i+1]])*points[[i+1]]$x[j]),
                 rad[[i+1]]*sin((1/rad[[i+1]])*points[[i+1]]$x[j]))
      
      vec1 <- ptsup - ptcentral
      vec2 <- ptinf - ptcentral
      
      v<-pmin(pmax(((vec1%*%vec2)[1,1])/
                     (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
      angles[j] <- acos(v)
    }
    tesener$Bending_energy[i] <- sum(alpha*((angles-pi)^2))/Lay
  }
  return(tesener)
}



gamad <- 0.15
lamad <- 0.04
alpha <- 1/100

cyl_width <- 5
cyl_length <- 20
xmin <- 0
ymin <- 0
xmax <- cyl_width
ymax <- cyl_length

rad_coef<-2.5

Radius <- cyl_width/(2*pi)
Radius2 <- rad_coef*Radius
cyl_width2 <- 2*pi*Radius2

cyl_thickness <- Radius2-Radius

n <- 100
Lay <- 20  #Layers
r <- cyl_width/n 
bet <- 10
steps <- 20

A0 <- (((Radius2+Radius)/2)*2*pi*cyl_length)/n

rec <- list()
rad <- list()
for(k in 1:Lay){
  rad[[k]]<- Radius+(k-1)*(cyl_thickness/(Lay-1)) #the radius of the layer k
  rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
}

ener<-energy_analysis(points)

p<-ggplot(ener, aes(x = Layer))+
  # geom_line(aes(y = Total_energy, colour = "Total energy"))+
  geom_line(aes(y=Tension_energy, colour = "Tension energy"))+
  geom_line(aes(y = Elastic_energy,colour = "Elastic energy"))+
  geom_line(aes(y = Contractile_energy, colour = "Contractile energy"))+
  geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
  labs(title = "Energy Analysis",
       x = "Layer", y = "Energy",
       color = "Energy type") +
  scale_colour_manual("",
                     breaks = c("Tension energy",
                                "Elastic energy",
                                "Contractile energy",
                                # "Total energy",
                                "Bending energy"),
                     values = c("red",
                                "blue",
                                # "darkgreen",
                                "orange",
                                "purple"))

show(p)