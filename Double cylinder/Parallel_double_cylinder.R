library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)


  move_points<-function(pt,wid,len,rc,n){
    ind<-sample(1:n,1)
    ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=rc)
    ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=rc)
    while((ptinx<0 || ptinx>wid)||(ptiny<0 || ptiny>len)){
      ptinx<-pt$x[ind]+rnorm(1,mean=0,sd=rc)
      ptiny<-pt$y[ind]+rnorm(1,mean=0,sd=rc)
    }
    pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+wid,ptinx+2*wid)
    pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
    return(pt)
  }
  
  tesellation_energy_double<-function(xt, yt, A0, rec1, rec2,
                                      R_c,gamad,lamad,n){
    tesel<-deldir(xt,yt,rw=rec1)
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perims<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas<-sapply(tilest,function(x){x$area})/A0
    tesener<-sum((areas-1)^2+(gamad/2)*(perims^2)+lamad*perims)
    
    tesel2<-deldir(xt*(R_c),yt,rw=rec2)
    tilest2<-tile.list(tesel2)[(n+1):(2*n)]
    perims2<-(tilePerim(tilest2)$perimeters)/sqrt(A0)
    areas2<-sapply(tilest2,function(x){x$area})/A0
    tesener2<-sum((areas2-1)^2+(gamad/2)*(perims2^2)+lamad*perims2)
    return(tesener+tesener2)
  }
  
  choice_metropolis<-function(delta,beta){
    p<-numeric()
    if(delta<=0){p<-1}else if(-delta*beta==0){p<-0}else{p<-exp(-delta*beta)}
    a=sample(c(0,1), size = 1, replace=TRUE, prob = c(1-p,p))
    return(a)
  }

  ggplotvor<-function(plotpoints,tit,wid){
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
  
  areasideplots<-function(px,py,rect,tit1,tit2,A0){
    tsl<-deldir(px,py,rw=rect)
    til<-tile.list(tsl)[(n+1):(2*n)]
    celledgearea<-data.frame(edges=numeric(),area=numeric())
    for (i in 1:length(til)) {
      celledgearea[i,c(1,2)]<-c(length(til[[i]]$x),til[[i]]$area/A0)
    }
    plotareaedges <- ggplot(celledgearea, aes(x = edges, y = area, colour = area))+
      geom_point()+xlab("Number of sides")+ylab("Relative area")+
      ggtitle(tit1)+
      stat_summary(aes(y = area,group=1), fun=mean, colour="#00BFC4", geom="line",group=1)
    show(plotareaedges)
    
    histedges <- ggplot(celledgearea,aes(edges))+geom_histogram(colour="#F8766D" ,fill="#7CAE00",bins=10,alpha=0.6)+
      xlab("Number of edges")+ylab("Frequency")+
      ggtitle(tit2)+
      xlim(2,11)
    show(histedges)
  }
  
  plotenergy <- function(en){
    ploten <- ggplot(en,aes(x=iteration,y=energy))+
      geom_line(colour="#F8766D")+
      xlab("Iteration of the algorithm")+
      ylab("Tesselation energy")+
      scale_y_continuous(trans = "log")+
      ggtitle("Energy relaxation of the tesselation")
    show(ploten)
  }
  
  #comienza el programa
  
  
  metropolisad<-function(seed = 666, steps = 250, n = 100,
                         RadiusA = 5/(2*pi), RadiusB = 2.5*5/(2*pi), cyl_length = 20,
                         gamma_ad = 0.15, lambda_ad = 0.04, beta = 100){
    
    
    #We define our variables
    
    cyl_width_A <- 2*pi*RadiusA
    cyl_width_B <- 2*pi*RadiusB
    
    #We define the vertices of the plane
    xmin <- 0
    xmax <- cyl_width_A
    ymin <- 0
    ymax <- cyl_length
    
    r <- cyl_width_A/n #radius to make the moves
    Am <- (cyl_width_A*(cyl_length))/n
    
    recA <- c(xmin,xmin+3*cyl_width_A,ymin,ymax)
    recB <- c(xmin,xmin+3*cyl_width_B,ymin,ymax)
    
    R_coef <- RadiusB/RadiusA
    
    
    #Start, first iteration
    
    set.seed(seed)
    
    x1 <- runif(n,xmin,xmax)
    y1 <- runif(n,ymin,ymax)
    x <-c(x1,x1+cyl_width_A,x1+2*cyl_width_A)
    y <-c(y1,y1,y1)
    
    points <- data.frame(x=x,y=y)
    pointsinit <- points
    energytesel <- tesellation_energy_double(points$x, points$y, Am, recA, recB,
                                             R_coef, gamma_ad, lambda_ad, n)
    energyinit <- energytesel
    
    #We create the variables to store the results
    energhist <- data.frame(iteration=numeric(steps), energy=numeric(steps))
    energhist[1,c(1,2)] <- c(0,energyinit)
    histpts <- data.frame(x=numeric(3*n*steps),
                          y=numeric(3*n*steps),
                          Frame = integer(3*n*steps))
    histpts[1:(3*n),c(1,2)] <- points
    histpts[1:(3*n),3] <- 1
    points2 <- data.frame(x=x,y=y)
    energytesel
    
    #Start of the loop
    for (j in 1:steps) {
      for(l in 1:n) {
        points2<-move_points(points,cyl_width_A,cyl_length,r,n)
        energytesel2<-tesellation_energy_double(points2$x, points2$y, Am,
                                                recA, recB, R_coef, 
                                                gamma_ad, lambda_ad, n)
        c<-choice_metropolis(energytesel2-energytesel,beta)
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
    return(list(histpts,energhist))
  }
  
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  
  results<-foreach(i=100:104, .combine = rbind, .packages = "deldir") %dopar% {
    metropolisad(seed = i, steps = 2)
  }
  
  stopImplicitCluster()
