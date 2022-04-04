library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(data.table)
library(stats)
library(nls.multstart)

changepoint3<-function(pt){
  ind<-sample(1:n,1)
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

energytes_adim<-function(tile){
  perim_ad<-(tilePerim(tile)$perimeters)/sqrt(A0)
  areas_ad<-sapply(tile,function(x){x$area})/A0
  energytesel<-sum((areas_ad-1)^2+(gam_ad/2)*perim_ad+lambda_ad*perim_ad)
  return(energytesel)
}

teselandenergy3_adim<-function(xt,yt){
  tesel<-deldir(xt,yt,rw=rec)
  tilest<-tile.list(tesel)[(n+1):(2*n)]
  perim_ad<-(tilePerim(tilest)$perimeters)/sqrt(A0)
  areas_ad<-sapply(tilest,function(x){x$area})/A0
  tesener<-sum((areas_ad-1)^2+(gam_ad/2)*perim_ad+lambda_ad*perim_ad)
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

byareaenergy<-function(points){
  tesel<-deldir(xt,yt,rw=rec)
  tilest<-tile.list(tesel)[(n+1):(2*n)]
  areas<-sapply(tilest,function(x){x$area})
  areas3<-c(areas,areas,areas)
  points[,3]<-areas3
  
  perim_ad<-(tilePerim(tilest)$perimeters)/sqrt(A0)
  areas_ad<-areas/A0
  encells<-(areas_ad-1)^2+(gam_ad/2)*perim_ad+lambda_ad*perim_ad
  encells3<-c(encells,encells,encells)
  
  names(points)<-c("x","y","area")
}

areasideplots<-function(xy){
  tsl<-deldir(xy$x,xy$y,rw=rec)
  til<-tile.list(tsl)[(n+1):(2*n)]
  celledgearea<-data.frame()
  for (i in 1:length(til)) {
    celledgearea[i,c(1,2)]<-c(length(til[[i]]$x),til[[i]]$area/A0)
  }
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

funaux<-function(p){
  tsl<-deldir(p$x,p$y,rw=rec)
  til<-tile.list(tsl)
  celledgearea<-data.frame(edges=integer(),area=double())
  for (i in 1:n) {
    celledgearea[i,c(1,2)]<-c(length(til[[i+n]]$x),til[[i+n]]$area/A0)
  }
  return (celledgearea)
}

funaux2<-function(histpt){
  ps<-250
  pts<-histpt[(350*3*n+1):((350+ps)*3*n),c(1,2,3)]
  edgear<-data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a<-length(edgear$edges)
    edgear[(a+1):(a+n),c(1,2)]<-funaux(pts[((i-1)*(3*n)+1):(i*3*n),c(1,2)])
    edgear[(a+1):(a+n),3]<-i
  }
  return(edgear)
}

stationarylewis<-function(edgear){
  #First execute funaux2 with the data, then this function makes the plots.
  plotareaedges<-ggplot(edgear, aes(x = edges, y = area, colour = area))+
    geom_point()+xlab("Number of sides")+ylab("Relative area")+
    ggtitle("Relative area of the cells by sides for the tessellations in stationary state")+
    stat_summary(aes(y = area,group=1), fun=mean, colour="#00BFC4", geom="line",group=1)
  show(plotareaedges)
  
  mined<-min(edgear$edges)
  maxed<-max(edgear$edges)
  quant<-maxed-mined
  
  datahist<-data.frame(edges=numeric(),frec=numeric())
  
  for (i in min(edgear$Frame):max(edgear$Frame)){
    for (j in mined:maxed) {
      dat<-filter(edgear, edges == j & Frame == i)
      pos<-(i-1)*quant+j-mined+1
      datahist[pos,c(1,2)]<-c(j, length(dat$edges))
    }
  }
  meandat<-data.frame(edges=numeric(), meanfrec=numeric())
  
  for (j in mined:maxed) {
    dat2<-filter(datahist, edges==j)
    pos<-j-mined+1
    meandat[pos,c(1,2)]<-c(j,mean(dat2$frec))
  }
  print(datahist)
  print(meandat)
  
  histedges<-ggplot(meandat,aes(edges,meanfrec))+
    geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
    xlab("Number of edges")+ylab("Average frequency")+
    ggtitle("Average quantity of sides of the cells for the tessellations in stationary state")+
    xlim(3,9)
  show(histedges)
}

regnls<-function(energh){
  x<-unlist(lapply(energh$iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls(y~I(a+b*(1-exp(-x/c))),start = list(a=2,b=0.3,c=0.1))
  plot(x,y)
  lines(x,predict(m),col="red",lwd=3)
  summary(m)
}
regnls2<-function(energh){
  x<-unlist(lapply(energh$iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls_multstart(y~I(a+b*(1-exp(-x/c))),
                   iter = 500,
                   start_lower = list(a=0,b=-5,c=0.01),
                   start_upper = list(a=5,b=5,c=1000))
  plot(x,y)
  lines(x,predict(m),col="red",lwd=3)
  summary(m)
  return(m)
}

plotenergy<-function(en){
  ploten<-ggplot(en,aes(x=iteration,y=energy))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Energy relaxation of the tesselation")
  show(ploten)
}

#comienza el programa

gam_ad<-0.15
lambda_ad<-0.04
xmin<-0
ymin<-0
xmax<-5
ymax<-20
wid<-xmax-xmin

n_adim<-100
n<-n_adim

r<- wid/n #radio en que cambiamos el punto
bet<-10
pasos<-250



A0<-(wid*(ymax-ymin))/n


x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x<-c(x1,x1+wid,x1+2*wid)
y<-c(y1,y1,y1)
points<-data.frame(x=x,y=y)

pointsinit<-points

rec<-c(xmin,xmin+3*wid,ymin,ymax)

teselacion <- deldir(points$x,points$y,rw=rec)
tiles <- tile.list(teselacion)
tilescyl<-tiles[(n+1):(2*n)]
energytesel<-energytes_adim(tilescyl)
energyinit<-energytesel
energhist<-data.frame(matrix(ncol=2,nrow=pasos))
names(energhist)<-c("iteration","energy")
energhist[1,c(1,2)]<-c(0,energyinit)
histpts<-data.frame(matrix(ncol=3,nrow=3*n*pasos))
names(histpts)<-c("x","y","Frame")
histpts[1:(3*n),c(1,2)]<-points
histpts[1:(3*n),3]<-1
points2<-data.frame(x=x,y=y)
energytesel

for (j in 1:pasos) {
  for(l in 1:n) {
    points2<-changepoint3(points)
    energytesel2<-teselandenergy3_adim(points2$x,points2$y)
    c<-choice(energytesel2-energytesel)
    cond<-c==1
    if(cond){
      points<-points2
      energytesel<-energytesel2
    }
  }
  gc()
  histpts[(j*3*n+1):(j*3*n+3*n),c(1,2)]<-points
  histpts[(j*3*n+1):(j*3*n+3*n),3]<-j+1
  energhist[j+1,c(1,2)]<-c(j,energytesel)
}
energyinit
energytesel


#plotvor3(points$x,points$y)

#plotvor3(pointsinit$x,pointsinit$y)

#ggplotvor(pointsinit, "       Initial Voronoi Tesselation")

#ggplotvor(points,"       Final Voronoi Tesselation")

#areasideplots(points)
energhist$energy<-energhist$energy/n #we plot mean cell energy
plotenergy(energhist)
