library(deldir)
library(ggplot2)
library(ggvoronoi)
library(dplyr)
library(stats)

move_points<-function(pt){
  ind<-sample(1:n,1)
  ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=r)
  ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=r)
  while((ptinx<xmin || ptinx>xmax)||(ptiny<ymin || ptiny>ymax)){
    ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=r)
    ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=r)
  }
  pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+cyl_width,ptinx+2*cyl_width)
  pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
  return(pt)
}


tesellation_energy_N<-function(xt,yt){
  tesener<-numeric(L)
  for (i in 1:L) {
    tesel<-deldir(xt*(rad[[i]]/Radius),yt,rw=rec[[i]])
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perims<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas<-sapply(tilest,function(x){x$area})/A0
    tesener[[i]]<-sum((areas-1)^2+(gamma_ad/2)*(perims^2)+lambda_ad*perims)
  }
  
  return(sum(tesener))
}

choice_metropolis<-function(delta){
  if(delta<=0){p<-1}else if(exp(-delta*bet)==Inf){p<-0}else{p<-exp(-delta*bet)}
  a=sample(c(0,1), size = 1,replace=TRUE, prob = c(1-p,p))
  return(a)
}

ggplotvor<-function(plotpoints,tit){
  rectangle <- data.frame(x=c(xmin,xmin,xmax+2*cyl_width,xmax+2*cyl_width),y=c(ymin,ymax,ymax,ymin))
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




gamma_ad <- 0.15
lambda_ad <- 0.04

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
L <- 10  #Layers

r <- cyl_width/n 
bet <- 10
steps <- 10


A0 <- (cyl_width*(ymax-ymin))/n


x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x <- c(x1,x1+cyl_width,x1+2*cyl_width)
y <- c(y1,y1,y1)
points <- data.frame(x=x,y=y)

pointsinit <- points

rec1 <- c(xmin,xmin+3*cyl_width,ymin,ymax)
rec2 <- c(xmin,xmin+3*cyl_width2,ymin,ymax)

rec <- list()
rad <- list()
for(k in 1:L){
  rad[[k]]<- Radius+(k-1)*(cyl_thickness/(L-1)) #the radius of the layer k
  rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
}

energytesel <- tesellation_energy_N(points$x,points$y)
energyinit <- energytesel
energhist <- data.frame(iteration = double(steps),energy = double(steps))
energhist[1,c(1,2)] <- c(0,energyinit)
histpts <- data.frame(x = double(3*n*steps), y = double(3*n*steps), Frame = double(3*n*steps))
histpts <- data.frame(matrix(ncol = 3,nrow = 3*n*steps))
names(histpts) <- c("x","y","Frame")
histpts[1:(3*n),c(1,2)] <- points
histpts[1:(3*n),3] <- 1
points2 <- data.frame(x = x, y = y)
energytesel

for (j in 1:steps) {
  for(l in 1:n) {
    points2 <- move_points(points)
    energytesel2 <- tesellation_energy_N(points2$x,points2$y)
    c <- choice_metropolis(energytesel2-energytesel)
    cond <- c==1
    if(cond){
      points <- points2
      energytesel <- energytesel2
    }
  }
  gc()
  histpts[(j*3*n+1):(j*3*n+3*n),c(1,2)] <- points
  histpts[(j*3*n+1):(j*3*n+3*n),3] <- j+1
  energhist[j+1,c(1,2)] <- c(j,energytesel)
}
save(list = c(histpts,energhist), file = "data300it.Rds")
energyinit
energytesel

energhist$energy <- energhist$energy/n
plot_energy(energhist)

#(points$x, points$y, rect = rec,
#              "Relative area of the cells by sides. Apical surface",
#              "Quantity of sides of the cells. Apical surface")
#areasideplots((Radius2/Radius)*points$x,points$y, rect = rec2,
#              "Relative area of the cells by sides. Basal surface",
#              "Quantity of sides of the cells. Basal surface")
