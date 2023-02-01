library(deldir)
library(ggplot2)
library(ggvoronoi)
library(dplyr)
library(stats)

# VERSION WITH BENDING  
#
# The main difference of this algorithm respect the original one is in the
# MOVEMENTS of the cells, since each layer moves separately.
# The other difference is in the ENERGY COMPUTATION, since we add a bending
# energy component, so the formula changes (and the function to compute it).
#
# All of this causes a change in the storage of the points.
# In the original algorithm, we just store the position of the apical layer, 
# since we know the position of each cell center in a layer by projections. In 
# this case, we have to store the position of the cell center in EVERY layer,
# because we don't project anymore. 
# That is why insted of a dataframe for the points, here we use a list of 
# dataframes.
# The rest of the code is similar to the original algorithm and works the same.
# Here we document the differences, to see the rest go to N_cylinder_algorithm.R
# in the N-cylinder folder.



bending_move_points<-function(pt){
  
  #We choose a cell
  
  ind<-sample(1:n,1)
  
  #Here we apply a movement to each layer. To do it efficiently we use lapply
  pt<- lapply(1:L, function(i){
    
    ptinx <- pt[[i]]$x[[ind]]+rnorm(1,mean=0,sd=r)
    ptiny <- pt[[i]]$y[[ind]]+rnorm(1,mean=0,sd=r)
    rd <- (rad[[i]]/Radius)
    # Check if the movement is inside the cylinder limits
    while((ptinx < xmin || ptinx > (rd*xmax)) ||
          (ptiny < ymin || ptiny > ymax)){
      ptinx <- pt[[i]]$x[[ind]]+rnorm(1, mean=0, sd=r)
      ptiny <- pt[[i]]$y[[ind]]+rnorm(1, mean=0, sd=r)
    }
    #Replicate the movement in the 2 copies of the cylinder rectangle
    pt[[i]]$x[c(ind,ind+n,ind+2*n)] <- c(ptinx, ptinx+rd*cyl_width, ptinx+2*rd*cyl_width)
    pt[[i]]$y[c(ind,ind+n,ind+2*n)] <- ptiny
    pt[[i]]
  })
  return(pt)
}

bending_tesellation_energy_N<-function(points){
  
  tesener<-numeric(L)
  
  #We compute energy as in the original algorithm
  
  tesener<-sapply(1:L,function(i) {
    tesel<-deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perims<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas<-sapply(tilest,function(x){x$area/A0})
    gam<-gamma_ad*exp((1-(rad[[i]]/rad[[1]]))/s0)
    sum((areas-1)^2+(gam/2)*(perims^2)+lambda_ad*perims)
  })
  
  #We add a bending energy component, depending on the angle that the center 
  # does with the center above and below itself. That is why the apical and 
  # basal centers don't have bending energy
  bendener <- sapply(2:(L-1), function(i){
    
    #First we compute the angles of every cell with the scalar product
    angles <- sapply(1:n, function(j){
      
      #we compute the position of a cell center in a layer of the 3d cylinder
      ptcentral <- c(points[[i]]$y[j],
                     rad[[i]]*cos((1/rad[[i]])*points[[i]]$x[j]),
                     rad[[i]]*sin((1/rad[[i]])*points[[i]]$x[j]))
      
      #we compute the position of the cell center bellow
      ptinf <- c(points[[i-1]]$y[j],
                 rad[[i-1]]*cos((1/rad[[i-1]])*points[[i-1]]$x[j]),
                 rad[[i-1]]*sin((1/rad[[i-1]])*points[[i-1]]$x[j]))
      
      #we compute the position of the cell center above
      ptsup <- c(points[[i+1]]$y[j],
                 rad[[i+1]]*cos((1/rad[[i+1]])*points[[i+1]]$x[j]),
                 rad[[i+1]]*sin((1/rad[[i+1]])*points[[i+1]]$x[j]))
      
      # we compute the scalar product and norm, with some adjustments so that 
      #the cos is between -1 and 1 (without pmin and pmax sometime we get a
      #slightly higher or lower result)
      vec1 <- ptsup - ptcentral
      vec2 <- ptinf - ptcentral
      v <- pmin(pmax(((vec1%*%vec2)[1,1])/
                       (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
      ang <- acos(v) #arcos to get angle
      return(ang) # we return the angles
    })
    
    #we use the angles to apply the bending energy formula
    return(sum(alpha*((angles-pi)^2)))
  })
  
  #we sum energies and return it, in this case the energy is normalized by layer.
  return((sum(tesener)+sum(bendener))/L)
}

choice_metropolis<-function(delta){
  if(delta<=0){p<-1}else if(exp(-delta*bet)==Inf){p<-0}else{p<-exp(-delta*bet)}
  a=sample(c(0,1), size = 1,replace=TRUE, prob = c(1-p,p))
  return(a)
}

ggplotvor<-function(plotpoints,tit){
  rectangle <- data.frame(x=c(xmin,xmin,xmax+2*cyl_width,xmax+2*cyl_width), y=c(ymin,ymax,ymax,ymin))
  pl <- ggplot(plotpoints,aes(x,y)) +
    geom_voronoi(aes(fill=as.factor(y)),size=.125, outline = rectangle,show.legend = FALSE) +
    geom_vline(xintercept = xmax,color = 'white',linetype='solid',size=1) +
    geom_vline(xintercept = xmax+cyl_width,color = 'white',linetype='solid',size=1) +
    stat_voronoi(geom="path",outline = rectangle) +
    geom_point(size=2) +
    theme(
      panel.grid.major = element_blank() # Remove gridlines (major)
      ,panel.grid.minor = element_blank() # Remove gridlines (minor)
      ,panel.background = element_blank() # Remove grey background
      ,plot.title = element_text(hjust = 0, size = 20, colour = "#323232") # Title size and colour
      ,plot.caption = element_text(vjust = 0.3, size = 11, colour = "#323232") # Caption size and colour
      ,axis.ticks.y = element_blank() # Remove tick marks (Y-Axis)
      ,axis.text.y =  element_blank() # Remove scale marks (Y-Axis)
      ,axis.title.y = element_blank() # Remove axis label (Y-Axis) 
      ,axis.ticks.x = element_blank() # Remove tick marks (X-Axis)
      ,axis.text.x  = element_blank() # Remove axis scale (X-Axis)
      ,axis.title.x = element_blank() # Remove axis label (X-Axis) 
      ,legend.position="bottom"
    ) +
    labs(title = tit # Title text
         ,caption = "Author: Eloy Serrano        ")
  show(pl)
}

areasideplots<-function(px, py, rect, tit, tit2){
  tsl<-deldir(px,py,rw=rect)
  til<-tile.list(tsl)[(n+1):(2*n)]
  celledgearea<-data.frame()
  for (i in 1:length(til)) {
    celledgearea[i,c(1,2)]<-c(length(til[[i]]$x),til[[i]]$area/A0)
  }
  colnames(celledgearea)<-c("edges","area")
  plotareaedges<-ggplot(celledgearea, aes(x = edges, y = area, colour = area))+
    geom_point()+xlab("Number of sides")+ylab("Relative area")+
    ggtitle(tit)+
    stat_summary(aes(y = area,group=1), fun=mean, colour="#00BFC4", geom="line",group=1)
  show(plotareaedges)
  
  histedges<-ggplot(celledgearea,aes(edges))+geom_histogram(colour="#F8766D", 
                                                            fill="#7CAE00",
                                                            bins=10,alpha=0.6)+
    xlab("Number of edges")+ylab("Frequency")+
    ggtitle(tit2)+
    xlim(2,11)
  show(histedges)
}

plot_energy<-function(en){
  ploten<-ggplot(en,aes(x=iteration,y=energy))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Energy relaxation of the tesselation")
  show(ploten)
}

nu_sq <- function(points, rec, n=100){
  Lay <- length(points)
  teselap <- deldir(points[[1]]$x, points[[1]]$y, rw = rec[[1]])
  teselba <- deldir(points[[Lay]]$x, points[[Lay]]$y, rw = rec[[Lay]])
  tilap <- tile.list(teselap)[(n+1):(2*n)]
  tilba <- tile.list(teselba)[(n+1):(2*n)]
  
  cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
  for (i in 1:length(tilap)) {
    cellsdf[i,c(1,2)]<-c(length(tilap[[i]]$x),length(tilba[[i]]$x))
  }
  num <- sum((cellsdf[,1]-cellsdf[,2])^2)/n
  den <- 2*(sum(cellsdf[,1])/n)*(sum(cellsdf[,2])/n)
  return(num/den)
}

gamma_ad <- 0.15
lambda_ad <- 0.04
alpha <- 1

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
L <- 3  #Layers

r <- cyl_width/n 
bet <- 10
steps <- 20

s0<-1

A0 <- ((Radius2+Radius)*pi*cyl_length)/n

rec <- list()
rad <- list()

for(k in 1:L){
  rad[[k]]<- Radius+(k-1)*(cyl_thickness/(L-1)) #the radius of the layer k
  rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
}

x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x <- c(x1,x1+cyl_width,x1+2*cyl_width)
y <- c(y1,y1,y1)

#We create a list of dataframes which will store the points in each layer.
points <- vector(mode = "list", length = L)
points <- lapply(1:L, function(i){data.frame(x=(rad[[i]]/Radius)*x,y=y)})

pointsinit <- points

energytesel <- bending_tesellation_energy_N(points)
energyinit <- energytesel
energhist <- data.frame(iteration=numeric(steps), energy=numeric(steps))
energhist[1,c(1,2)] <- c(0,energyinit)


# Also the storage of the historic and evolution of the algorithm changes, 
# now we have a list of dataframes, each element of the list is, in order,
# an iteration of the algorithm (1st element correspond to the first iteration,
# 2nd element to the 2nd iteration, and so on)
histpts <- vector(mode="list", length = L)

for (i in 1:L) {
  histpts[[i]] <- data.frame(x = double(3*n*steps), y = double(3*n*steps), Frame = double(3*n*steps))
}

for (i in 1:L) {
  histpts[[i]][1:(3*n),c(1,2)] <- points[[i]]
  histpts[[i]][1:(3*n),3] <- i+1
}



# ALGORITHM STARTS


for (j in 1:steps) {
  for(l in 1:n) {
    points2 <- bending_move_points(points)
    energytesel2 <- bending_tesellation_energy_N(points2)
    c <- choice_metropolis(energytesel2-energytesel)
    cond <- c==1
    if(cond){
      points <- points2
      energytesel <- energytesel2
    }
    rm(points2, energytesel2)
    gc()
  }
  
  #We store the points of each layer after each iteration
  
  for (i in 1:L) {
    histpts[[i]][(j*3*n+1):(j*3*n+3*n),c(1,2)] <- points[[i]]
    histpts[[i]][(j*3*n+1):(j*3*n+3*n),3] <- j+1
  }
  
  #We store the energy of the tessellation after each iteration
  
  energhist[j+1,c(1,2)] <- c(j,energytesel)
  
}

save(histpts,energhist, file = "data300it.Rds")

nu2 <- nu_sq(points = points, rec = rec, n = n)

energyinit
energytesel

energhist$energy <- energhist$energy/n
plot_energy(energhist5)

#(points$x, points$y, rect = rec,
#              "Relative area of the cells by sides. Apical surface",
#              "Quantity of sides of the cells. Apical surface")
#areasideplots((Radius2/Radius)*points$x,points$y, rect = rec2,
#              "Relative area of the cells by sides. Basal surface",
#              "Quantity of sides of the cells. Basal surface")
