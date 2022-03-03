library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)


changepoint3_alg2<-function(pt,ind){
  ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=r)
  ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=r)
  while((ptinx<xmin || ptinx>xmax)||(ptiny<ymin || ptiny>ymax)){
    ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=r)
    ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=r)
  }
  pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+wid,ptinx+2*wid)
  pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
  return(pt)
}

energytes<-function(tile){
  perim<-tilePerim(tile)$perimeters
  energycell<-c()
  for (i in 1:length(tile)){
    energycell[i]<-((k/2)*(((tile[[i]]$area)-A0)^2)) #+(gam/2)*(perim[[i]])^2+lambda*perim[[i]]
  }
  energytesel<-sum(energycell)
  return(energytesel)
}

teselandenergy3<-function(xt,yt){
  tesel<-deldir(xt,yt,rw=rec)
  tilest<-tile.list(tesel)[(n+1):(2*n)]
  tesener<-energytes(tilest)
  return(tesener)
}

choice<-function(delta){
  if(delta<=0){p<-1}else if(exp(-delta*bet)==Inf){p<-0}else{p<-exp(-delta*bet)}
  a=sample(c(0,1), size = 1,replace=TRUE, prob = c(1-p,p))
  return(a)
}

plotvor3<-function(xv,yv){
  tesv<-deldir(xv,yv,rw=rec)
  tilesv<-tile.list(tesv)
  plot(tilesv, pch=19, border="white", fillcol=hcl.colors(10,"Blues"))
  lines(c(xmax,xmax),c(ymin,ymax))
  lines(c(xmax+wid,xmax+wid),c(ymin,ymax))
  plot(tilesv[(n+1):(2*n)], pch=19, border="white", fillcol=hcl.colors(10,"Blues"))
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

plotminener<-function(histpt){
  #Lets draw the tesselation corresponding with minimum energy
  energhist<-list()
  for (i in 1:length(histpt)) {
    energhist[[i]]<-histpt[[i]][[2]]
  }
  indmin<-which.min(unlist(energhist))
  ptsmin<-histpt[[indmin]][[1]]
  plotvor3(ptsmin$x,ptsmin$y)
}

trhistpts<-function(histpt){
  #prepara los puntos para hacer el video
  t<-length(histpt)
  df1<-histpt[[1]][[1]]
  df1[,3]<-1
  names(df1)<-c("x","y","Frame")
  for (i in 2:t) {
    a<-length(df1$x)
    df1[a+1:(3*n),c(1,2)]<-histpt[[i]][[1]]
    df1[a+1:(3*n),3] <- i
  }
  return(df1)
}

#comienza el programa

k<-1
gam<-1
lambda<-1
xmin<-0
ymin<-0
xmax<-10
ymax<-10
wid<-xmax-xmin
n<-35

r<- wid/n #radio en que cambiamos el punto
bet<-10
pasos<-50


A0<-(wid*(ymax-ymin))/n


x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x<-c(x1,x1+wid,x1+2*wid)
y<-c(y1,y1,y1)
points<-data.frame(x=x,y=y)

pointsinit<-points

#histpts<-rep(list(list()), pasos)
#histpts<-vector(mode="logical",length=pasos)
#histpts[[1]]<-c(x,y)
#para guardar los puntos que se van generando

rec<-c(xmin,xmin+3*wid,ymin,ymax)

teselacion <- deldir(points$x,points$y,rw=rec)
tiles <- tile.list(teselacion)
tilescyl<-tiles[(n+1):(2*n)]
energytesel<-energytes(tilescyl)
energyinit<-energytesel
energhist<-data.frame(iteration=0,energy=energyinit)

histpts<-list()

#añadir un plot aqui

for (j in 1:pasos) {
  for(l in 1:n) {
    points2<-changepoint3_alg2(points,l)
    energytesel2<-teselandenergy3(points2$x,points2$y)
    c<-choice(energytesel2-energytesel)
    if(c==1){
      points<-points2
      energytesel<-energytesel2
      histpts[[j*l]]<-list(points,energytesel)
    }
  }
  energhist[j+1,c(1,2)]<-c(j,energytesel)
}
energyinit
energytesel
histpts<-Filter(Negate(is.null), histpts)


#plotvor3(points$x,points$y)

#plotvor3(pointsinit$x,pointsinit$y)

ggplotvor(pointsinit, "       Initial Voronoi Tesselation")

ggplotvor(points,"       Final Voronoi Tesselation")

areasideplots(points)
plotenergy(energhist)