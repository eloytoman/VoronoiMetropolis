library(deldir)
library(ggplot2)
library(ggvoronoi)
library(stats)
library(dplyr)
library(nls.multstart)


#This script takes the results after performing dinstinct simulations
#and makes the plots corresponding to the lewis law and the histogram of 
#sides of al the 400ths tessellations of every simulation

funaux<-function(p){
  tsl<-deldir(p$x,p$y,rw=rec)
  til<-tile.list(tsl)
  celledgearea<-data.frame(edges=integer(),area=double())
  for (i in 1:n) {
    celledgearea[i,c(1,2)]<-c(length(til[[i+n]]$x),til[[i+n]]$area/A0)
  }
  return (celledgearea)
}

funaux2sim<-function(ptsord){
  #ps is the number of simulations done
  ps<-1000
  edgear<-data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a<-length(edgear$edges)
    edgear[(a+1):(a+n),c(1,2)]<-funaux(ptsord[((i-1)*(3*n)+1):(i*3*n),c(1,2)])
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
    xlim(2,11)
  show(histedges)
}

ord<-function(results){
  #With this function we extract the points for the 400th iteration of every simulation
  
  ptsord<-data.frame(x=numeric(300*1000),y=numeric(300*1000))
  for (i in 0:999) {
    ptsord[(i*300+1):((i+1)*300),c(1,2)]<-results[[i+1]][results[[i+1]]$Frame==300,c(1,2)]
  }
  return(ptsord)
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
  return(summary(m)[["coefficients"]][c(1,2,3)])
}

adjsim<-function(results){
  #With this function we extract the points for the 300th iteration of every simulation
  coefest<-data.frame(a=double(1000),b=double(1000),c=double(1000))
  for (i in 1:1000) {
    coefest[i,c(1,2,3)]<-regnls2(results[[i+1001]])
  }
  a<-mean(coefest$a)
  b<-mean(coefest$b)
  c<-mean(coefest$c)
  adj<-function(x) a+b*(1-exp(-x/c))
  adj_data<-data.frame(x=1:300,y=adj(1:300))
  ploten<-ggplot(data = adj_data, aes(x=x,y=y))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Average energy relaxation of the tesselation. 1000 simulations")
  show(ploten)
  return(c(a,b,c))
}


scutoids_analysis_oneiter<-function(pointsAx, pointsAy, pointsBx, pointsBy, rect1, rect2, n = 100){
  tslA<-deldir(pointsAx,pointsAy,rw=rect1)
  tilA<-tile.list(tslA)[(n+1):(2*n)]
  tslB<-deldir(pointsBx,pointsBy,rw=rect2)
  tilB<-tile.list(tslB)[(n+1):(2*n)]
  
  cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
  
  for (i in 1:length(tilA)) {
    cellsdf[i,c(1,2)]<-c(length(tilA[[i]]$x),length(tilB[[i]]$x))
  }
  countdf<- cellsdf %>%
    group_by(edgesA,edgesB) %>%
    summarize(count=n())
  
  scutoidsplot<-ggplot(countdf, aes(x = edgesA, y = edgesB, label=count))+
    geom_count(shape = "square", aes(color= count))+
    xlab("Edges on apical surface")+ylab("Edges on basal surface")+
    ggtitle("Polygon class of apical and basal surfaces")+
    guides(colour = "colorbar", size = "none")+
    scale_size_area(max_size = 30)+
    geom_label()+
    scale_fill_gradient(low = "light blue", high = "deepskyblue")+
    xlim(3.3,7.5)+
    ylim(3.5,7.5)
  show(scutoidsplot)
  
  return(countdf)
}



scutoids_prep <- function(pointsAx,pointsAy,pointsBx,pointsBy,rect1,rect2,n){
  tslA <- deldir(pointsAx,pointsAy,rw=rect1)
  tilA <- tile.list(tslA)[(n+1):(2*n)]
  tslB <- deldir(pointsBx,pointsBy,rw=rect2)
  tilB <- tile.list(tslB)[(n+1):(2*n)]
  cellsdf <- data.frame(edgesA=integer(),edgesB=integer())
  for (i in 1:length(tilA)) {
    cellsdf[i,c(1,2)] <- c(length(tilA[[i]]$x),length(tilB[[i]]$x))
  }
  countdf <- cellsdf %>%
    group_by(edgesA,edgesB) %>%
    summarize(count=n())
  return(countdf)
}


scutoids_analysis_stationary <- function(histpts, rect1, rect2, n = 100){
  lon <- 50 #how many iterations we want to have
  histdf_count <- data.frame(edgesA=integer(lon*n),edgesB=integer(lon*n),count=integer(lon*n))
  for (i in 1:lon) {
    histdf_count[((i-1)*n+1):(i*n)]<-
      scutoids_prep(filter(histpts, Frame==250+i-1)$x,
                    filter(histpts, Frame==250+i-1)$y,
                    (Radius2/Radius)*filter(histpts, Frame==250+i-1)$x,
                    filter(histpts, Frame==250+i-1)$y,
                    rect1,rect2)
  }
  histdf_avgcount <- histdf_count %>%
    group_by(edgesA,edgesB) %>%
    summarize(avg_count=mean(count))
  
  
  scutoidsplot<-ggplot(histdf_avgcount, aes(x = edgesA, y = edgesB, label=avg_count))+
    geom_count(shape = "square", aes(color= avg_count))+
    xlab("Average of edges on apical surface")+ylab("Average of edges on basal surface")+
    ggtitle("Average polygon class of apical and basal surfaces")+
    guides(colour = "colorbar", size = "none")+
    scale_size_area(max_size = 30)+
    geom_label()+
    scale_fill_gradient(low = "light blue", high = "deepskyblue")+
    xlim(3.3,7.5)+
    ylim(3.5,7.5)
  show(scutoidsplot)
  return(histdf_avgcount)
}

resord<-ord(results1000)
edgearsim<-funaux2sim(resord)
stationarylewis(edgearsim[1:1000,c(1,2,3)])
coef<-adjsim(results1000)

scutoids_analysis_oneiter(points$x,points$y,(Radius2/Radius)*points$x,points$y,rec,rec2)



