library(deldir)
library(ggvoronoi)
library(ggplot2)
library(plotly)
library(av)



energytes<-function(tile){
  perim<-tilePerim(tile)$perimeters
  energytesel<-c(0)
  for (i in 1:length(tile)){
    energytesel<-energytesel+(k/2)*(tile[[i]]$area-A0)^2+
      (gam/2)*(perim[[i]])^2+lambda*perim[[i]]
  }
  return(energytesel)
}

teselandenergy3<-function(xt,yt){
  tesel<-deldir(xt,yt,rw=rec)
  tilest<-tile.list(tesel)[(n+1):(2*n)]
  perim<-tilePerim(tilest)$perimeters
  tesener<-c(0)
  for (i in 1:n){
    tesener<-tesener+(k/2)*(tilest[[i]]$area-A0)^2+
      (gam/2)*(perim[[i]])^2+lambda*perim[[i]]
  }
  
  return(tesener)
}

plotvor3<-function(xv,yv){
  tesv<-deldir(xv,yv,rw=rec)
  tilesv<-tile.list(tesv)
  plot(tilesv, pch=19, border="white", fillcol=hcl.colors(10,"Greens"))
  lines(c(xmax,xmax),c(ymin,ymax))
  lines(c(xmax+wid,xmax+wid),c(ymin,ymax))
  plot(tilesv[(n+1):(2*n)], pch=19, border="white", fillcol=hcl.colors(10,"Greens"))
  lines(c(xmax,xmax),c(ymin,ymax))
  lines(c(xmax+wid,xmax+wid),c(ymin,ymax))
}

ggplotvor<-function(plotpoints,tit){
  rectangle <- data.frame(x=c(xmin,xmin,xmax+2*wid,xmax+2*wid),y=c(ymin,ymax,ymax,ymin))
  pl <- ggplot(plotpoints,aes(x,y)) +
    geom_voronoi(aes(fill=as.factor(y)),size=.125, outline = rectangle,show.legend = FALSE) +
    geom_vline(xintercept = xmax,color = 'white',linetype='solid',size=1) +
    geom_vline(xintercept = xmax+wid,color = 'white',linetype='solid',size=1) +
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

areasideplots<-function(xy){
  tsl<-deldir(xy$x,xy$y,rw=rec)
  til<-tile.list(tsl)
  #Lwlw<-lawSummary(tsl)
  #option1
  celledgearea<-data.frame()
  for (i in 1:length(til)) {
    celledgearea[i,c(1,2)]<-c(length(til[[i]]$x),til[[i]]$area/A0)
  }
  #option2
  #celledgearea<-data.frame(Lwlw$num.edges,Lwlw$tile.areas/A0)
  colnames(celledgearea)<-c("edges","area")
  plotareaedges<-ggplot(celledgearea, aes(x = edges, y = area, colour = area))+
    geom_point()+xlab("Number of sides")+ylab("Relative area")+
    ggtitle("Relative area of the cells by sides")+
    stat_summary(aes(y = area,group=1), fun=mean, colour="#00BFC4", geom="line",group=1)
  show(plotareaedges)
  
  histedges<-ggplot(celledgearea,aes(edges))+geom_histogram(colour="#F8766D" ,fill="#7CAE00",bins=10,alpha=0.6)+
    xlab("Number of edges")+ylab("Frequency")+
    ggtitle("Quantity of sides of the cells")+
    xlim(2,11)
  show(histedges)
  
}

plotenergy<-function(en){
  ploten<-ggplot(en,aes(x=iteration,y=energy))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Tesselation energy")+
    ggtitle("Energy relaxation of the tesselation")
  show(ploten)
}

lloydpoints<-function(x){
  if(x-wid>xmax){x<-x-wid}else if(x-wid<xmin){x<-x+wid}else{x}
}

plotcompener<-function(en1,en2,name1,name2){
  #After introducing the energy vector of lloyd and metropolis algorithm,
  #this function makes the combined plot of both energies.
  m=length(en1$iteration)
  energy<-en1
  energy[(m+1):(2*m),c(1,2)]<-en2
  energy[1:m,3]<-name1
  energy[(m+1):(2*m),3]<-name2
  colnames(energy)<-c("iteration","energy","algorithm")
  
  p<-ggplot(energy,aes(x=iteration, y=energy, color=algorithm))+
    geom_line()+
    xlab("Iteration of the algorithm")+
    ylab("Tesselation energy")+
    ggtitle("Energy relaxation of the tesselation by algorithms")
  show(p)
}


#comienza el programa

k<-1
gam<-1
lambda<-1
xmin<-0
ymin<-0
xmax<-1
ymax<-1
wid<-xmax-xmin
n<-35

r<- wid*2/n #radio en que cambiamos el punto
bet<-1000
pasos<-50


A0<-(wid*(ymax-ymin))/n


x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x<-c(x1,x1+wid,x1+2*wid)
y<-c(y1,y1,y1)
points<-data.frame(x=x,y=y)

ptsinitx<-x
ptsinity<-y
pointsinit<-data.frame(x=ptsinitx,y=ptsinity)


rec<-c(xmin,xmin+3*wid,ymin,ymax)

teselacion <- deldir(points$x,points$y,rw=rec)
tiles <- tile.list(teselacion)
tilescyl<-tiles[(n+1):(2*n)]
energytesel<-energytes(tilescyl)
energyinit<-energytesel
energhist<-data.frame(iteration=0,energy=energyinit)

for (j in 1:pasos) {
  teselacion <- deldir(points$x,points$y,rw=rec)
  tiles <- tile.list(teselacion)
  tilescyl<-tiles[(n+1):(2*n)]
  energytesel<-energytes(tilescyl)
  energhist[j+1,c(1,2)]<-c(j,energytesel)
  
  centroids<-tile.centroids(tiles)[(n+1):(2*n),c(1,2)]
  centroids$x<-unlist(sapply(centroids$x,lloydpoints))
  centroids1<-centroids
  centroids1$x<-centroids1$x-wid
  centroids2<-centroids
  centroids2$x<-centroids$x+wid
  points[1:n,c(1,2)]<-centroids1
  points[(n+1):(2*n),c(1,2)]<-centroids
  points[(2*n+1):(3*n),c(1,2)]<-centroids2
}

energyinit
energytesel
#histpts<-Filter(Negate(is.null), histpts)

ggplotvor(pointsinit,"       Initial Voronoi Tesselation (Lloyd)")

ggplotvor(points,"       Final Voronoi Tesselation (Lloyd)")

areasideplots(points)
plotenergy(energhist)





