library(deldir)
library(ggplot2)
library(ggvoronoi)
library(dplyr)
library(stats)


#FUNCTIONS INVOLVED IN THE ALGORITHM

#moving points function

move_points<-function(pt){
  ind<-sample(1:n,1) #we choose randomly one of the cells
  #Now we do a random movement of the center of the cell, following a gaussian
  #distribution with mean 0 and std = r (r param is defined below)
  ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=r)
  ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=r)
  
  #Now we check that the movement doesn't push the center outside the rectangle
  while((ptinx<xmin || ptinx>xmax)||(ptiny<ymin || ptiny>ymax)){
    ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=r)
    ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=r)
  }
  
  #If the movement is correct, we replicate it (remember that we have 3 
  # copies of the rectangle, to better simulate the cylinder's unrolling edge)
  pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+cyl_width,ptinx+2*cyl_width)
  pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
  return(pt)
}


#Function to compute the energy of the tessellation

tesellation_energy_N<-function(xt,yt){
  
  #Intialite the vector of the energy, L components for each of the L layers
  tesener<-numeric(L)
  
  #We use sapply to compute in a vectorized way the energy in the L layers 
  #This is more efficient since for algorithms take a lot of computing resources 
  #and time
  
  tesener<-sapply(1:L,function(i) {
    
    #deldir computes the voronoi tessellation of all the cells
    #   xt is the vector with all the horizontal coordinates of the cells
    #   yt is the vector with the vertical coordinates of the cells
    #   rw is the rectangle which is the boundary for the tessellation.
    #   we multiply by a constant to obtain the exact coordinate in each layer
    tesel<-deldir(xt*(rad[[i]]/Radius),yt,rw=rec[[i]])
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perims<-(tilePerim(tilest)$perimeters)/sqrt(A0) #we compute the perimeters
    areas<-sapply(tilest,function(x){x$area})/A0 #we compute the areas
    
    #we apply the energy formula
    sum((areas-1)^2+(gamma_ad/2)*(perims^2)+lambda_ad*perims)
  })
  #we return the total energy of the tessellation
  return(sum(tesener)/L)
}


#Function to compare energies

choice_metropolis<-function(delta){
  #delta is the energy difference, if it is negative, we make the move with
  #probability 1, if it is positive, we do it with the specified probability
  # the case where the exponential is infinite is to avoid some errors
  if(delta<=0){p<-1}else if(exp(-delta*bet)==Inf){p<-0}else{p<-exp(-delta*bet)}
  
  #we use the probabiliy to obtain 1 if we do the move, 0 if not
  a=sample(c(0,1), size = 1,replace=TRUE, prob = c(1-p,p))
  return(a)
}


#function to plot the voronoi cells

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


#function to make a histogram of cell edges and its frequency
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


#Function to plot the energy relaxation

plot_energy<-function(en){
  ploten<-ggplot(en,aes(x=iteration,y=energy))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Energy relaxation of the algorithm")
  show(ploten)
}


#function that computes n^2 parameter
nu_sq <- function(points, rec, RadB= 2.5*5/(2*pi) , n=100){
  Lay <- length(rec)
  teselap <- deldir(points$x, points$y, rw = rec[[1]])
  teselba <- deldir(RadB*points$x, points$y, rw = rec[[Lay]])
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



# PARAMETERS OF THE ALGORITHM

#Energy parameters

gamma_ad <- 0.15
lambda_ad <- 0.04


#Cylinder dimensions

cyl_width <- 5 #width of the rectangle before unrolling
cyl_length <- 20 #length of cylinder
xmin <- 0
ymin <- 0
xmax <- cyl_width
ymax <- cyl_length

rad_coef<-2.5

Radius <- cyl_width/(2*pi) #Apical radius
Radius2 <- rad_coef*Radius #Basal radius
cyl_width2 <- 2*pi*Radius2

cyl_thickness <- Radius2-Radius

n <- 100 #Cells
L <- 10  #Layers


#Metropolis algorithm parameters

r <- cyl_width/n #movement parameter
bet <- 10  #acceptance parameter
steps <- 50  #steps of the algorithm



A0 <- ((Radius2+Radius)*pi*cyl_length)/n  #Ideal area


#GENERATION OF THE INITIAL CELL CONFIGURATION

x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x <- c(x1,x1+cyl_width,x1+2*cyl_width)
y <- c(y1,y1,y1)
points <- data.frame(x=x,y=y)

pointsinit <- points

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
histpts[1:(3*n),c(1,2)] <- points
histpts[1:(3*n),3] <- 1
points2 <- data.frame(x = x, y = y)
energytesel


# START OF THE ALGORITM

for (j in 1:steps) {
  for(l in 1:n) {
    
    #First we do a movement
    points2 <- move_points(points)
    
    #We compute its energy
    energytesel2 <- tesellation_energy_N(points2$x,points2$y)
    
    # We compare to previous energy
    c <- choice_metropolis(energytesel2-energytesel)
    
    # If the energy is lower, we accept the movement, if it is not, we accept
    # the movement with certain probability
    
    cond <- c==1
    if(cond){
      points <- points2
      energytesel <- energytesel2
    }
  }
  
  gc() # Memory clear
  
  #We save the results of the iteration to have the historic
  histpts[(j*3*n+1):(j*3*n+3*n),c(1,2)] <- points
  histpts[(j*3*n+1):(j*3*n+3*n),3] <- j+1
  energhist[j+1,c(1,2)] <- c(j,energytesel)
}

#We save the results
save(histpts,energhist, file = "data300it.Rds")

#We compute the nu^2 parameter to obtain more information of the final topology
nu2 <- nu_sq(points = points, rec = rec, n = n)


energyinit
energytesel

#normalization of the energy by cell and plot
energhist$energy <- energhist$energy/n
plot_energy(energhist)


#Plot of the topology of the cell

#(points$x, points$y, rect = rec,
#              "Relative area of the cells by sides. Apical surface",
#              "Quantity of sides of the cells. Apical surface")
#areasideplots((Radius2/Radius)*points$x,points$y, rect = rec2,
#              "Relative area of the cells by sides. Basal surface",
#              "Quantity of sides of the cells. Basal surface")
