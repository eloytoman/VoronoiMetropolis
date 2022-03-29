library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)


system.time({
  changepoint3<-function(pt,xmin,xmax,ymin,ymax,rc,n){
    wid<-xmax-xmin
    ind<-sample(1:n,1)
    ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=rc)
    ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=rc)
    while((ptinx<xmin || ptinx>xmax)||(ptiny<ymin || ptiny>ymax)){
      ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=rc)
      ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=rc)
    }
    pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+wid,ptinx+2*wid)
    pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
    return(pt)
  }
  
  energytes_adim<-function(tile,A0,gamad,lamad){
    perim_ad<-(tilePerim(tile)$perimeters)/sqrt(A0)
    areas_ad<-sapply(tile,function(x){x$area})/A0
    energytesel<-sum((areas_ad-1)^2+(gamad/2)*perim_ad+lamad*perim_ad)
    return(energytesel)
  }
  
  teselandenergy3_adim<-function(xt,yt,A0,rect,gamad,lamad,n){
    tesel<-deldir(xt,yt,rw=rect)
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perim_ad<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas_ad<-sapply(tilest,function(x){x$area})/A0
    tesener<-sum((areas_ad-1)^2+(gamad/2)*perim_ad+lamad*perim_ad)
    return(tesener)
  }
  
  choice<-function(delta,beta){
    if(delta<=0){p<-1}else if(exp(-delta*beta)==Inf){p<-0}else{p<-exp(-delta*beta)}
    a=sample(c(0,1), size = 1,replace=TRUE, prob = c(1-p,p))
    return(a)
  }
  
  byareaenergy<-function(points,rec,n,A0,gam_ad,lambda_ad){
    tesel<-deldir(xt,yt,rw=rec)
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    areas<-sapply(tilest,function(x){x$area})
    areas3<-c(areas,areas,areas)
    points[,3]<-areas3
    
    perim_ad<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas_ad<-sapply(tilest,function(x){x$area})/A0
    encells<-(areas_ad-1)^2+(gam_ad/2)*perim_ad+lambda_ad*perim_ad
    encells3<-c(encells,encells,encells)
    
    names(points)<-c("x","y","area")
  }
  
  plotvor3<-function(xv,yv,rec){
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
  
  areasideplots<-function(xy,rec,A0){
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
      scale_y_continuous(trans = "log")+
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
      scale_y_continuous(name="Tesselation energy", trans="log")+
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
  
  #comienza el programa
  
  
  metropolisad<-function(seed=666,pasos=250,n_adim=100,
                         xmin=0,xmax=5,ymin=0,ymax=20,
                         gam_ad,lambda_ad,bet){
    wid<- xmax-xmin
    n <- n_adim
    r <- wid/n #radio en que cambiamos el punto
    Am <- (wid*(ymax-ymin))/n
    set.seed(seed)
    x1 <- runif(n,xmin,xmax)
    y1 <- runif(n,ymin,ymax)
    x <-c(x1,x1+wid,x1+2*wid)
    y <-c(y1,y1,y1)
    points <- data.frame(x=x,y=y)
    pointsinit <- points
    
    rec <- c(xmin,xmin+3*wid,ymin,ymax)
    
    teselacion <- deldir(points$x,points$y,rw=rec)
    tiles <- tile.list(teselacion)
    tilescyl <- tiles[(n+1):(2*n)]
    
    energytesel <- energytes_adim(tilescyl,Am,gam_ad,lambda_ad)
    energyinit <- energytesel
    energhist <- data.frame(matrix(ncol=2,nrow=pasos))
    names(energhist) <- c("iteration","energy")
    energhist[1,c(1,2)] <- c(0,energyinit)
    histpts <- data.frame(matrix(ncol=3,nrow=3*n*pasos))
    names(histpts) <- c("x","y","Frame")
    histpts[1:(3*n),c(1,2)] <- points
    histpts[1:(3*n),3] <- 1
    points2<-data.frame(x=x,y=y)
    energytesel
    
    
    for (j in 1:pasos) {
      for(l in 1:n) {
        points2<-changepoint3(points,xmin,xmax,ymin,ymax,rec,n_adim)
        energytesel2<-teselandenergy3_adim(points2$x,points2$y,Am,rec,
                                           gam_ad,lambda_ad,n)
        c<-choice(energytesel2-energytesel,bet)
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
      energytesel
    }
    #plotvor3(points$x,points$y)
    
    #plotvor3(pointsinit$x,pointsinit$y)
    
    #ggplotvor(pointsinit, "       Initial Voronoi Tesselation")
    
    #ggplotvor(points,"       Final Voronoi Tesselation")
    
    #areasideplots(points,Am)
    #we plot mean cell energy
    return(list(histpts,energhist))
  }
  
  #metropolisad(seed = 123, pasos = 5, n_adim = 100,
  #             xmin = 0, xmax = 5, ymin = 0, ymax = 20, 
  #             gam_ad = 0.15, lambda_ad = 0.04, bet=10)
  
  
  
  #cl <- makeCluster(4)
  #registerDoParallel(cl)
  
  
  
  #cl <- makeCluster(4)
  #registerDoParallel(cl)
  
  
  for(i in 100:104) {
    metropolisad(seed = i, pasos = 20, n_adim = 100,
                 xmin = 0, xmax = 5, ymin = 0, ymax = 20,
                 gam_ad = 0.15, lambda_ad = 0.04, bet = 10)
  }
 # results<-foreach(i=100:104, .combine = rbind, .packages = "deldir") %dopar% {
  #  metropolisad(seed = i, pasos = 20, n_adim = 100,
   #              xmin = 0, xmax = 5, ymin = 0, ymax = 20,
    #             gam_ad = 0.15, lambda_ad = 0.04, bet = 10)
  #}
  
  #stopImplicitCluster()
})