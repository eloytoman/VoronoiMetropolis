library(deldir)
library(ggplot2)
#library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)


  bending_move_points<-function(pt, wid, len, rc, n = 100, Lay = 3, rad){
    
    ind <-sample(1:n,1)
    
    pt <- lapply(1:Lay, function(i){
      
      ptinx <- pt[[i]]$x[[ind]]+rnorm(1,mean=0,sd=rc)
      ptiny <- pt[[i]]$y[[ind]]+rnorm(1,mean=0,sd=rc)
      s <- (rad[[i]]/rad[[1]])
      while((ptinx < 0 || ptinx > (s*wid)) ||
            (ptiny < 0 || ptiny > len)){
        ptinx <- pt[[i]]$x[[ind]] + rnorm(1, mean=0, sd=rc)
        ptiny <- pt[[i]]$y[[ind]] + rnorm(1, mean=0, sd=rc)
      }
      pt[[i]]$x[c(ind,ind+n,ind+2*n)] <- c(ptinx, ptinx+s*wid, ptinx+2*s*wid)
      pt[[i]]$y[c(ind,ind+n,ind+2*n)] <- ptiny
      pt[[i]]
    })
    return(pt)
  }
  
  bending_tesellation_energy_N <- function(points, A0, rec, rad, gamad, lamad,
                                           alpha = 1, n, Lay = 3, s0=1){
    
    tesener <- numeric(L)
    
    tesener <- sapply(1:Lay, function(i){
      
      tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
      tilest <- tile.list(tesel)[(n+1):(2*n)]
      
      perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
      areas <- sapply(tilest,function(x){x$area/A0})
      
      gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/s0)
      
      sum((areas-1)^2+(gam/2)*(perims^2)+lamad*perims)
    })
    
    bendener <- sapply(2:(Lay-1), function(i){
      
      angles <- sapply(1:n, function(j){
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
        
        #We use pmin and pmax to avoid errors in the arc-cosine computation
        v <- pmin(pmax(((vec1%*%vec2)[1,1])/
                       (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
        ang <- acos(v)        
        return(ang)
      })
      
      return(sum(alpha*((angles-pi)^2)))
    })
    
    return((sum(tesener)+sum(bendener))/Lay)
  }
  
  choice_metropolis<-function(delta,beta){
    
    p<-numeric(1)
    
    if(delta<=0){p<-1}else if(exp(-delta*beta)==Inf){p<-0}else{p<-exp(-delta*beta)}
    
    a=sample(c(0,1), size = 1, replace=TRUE, prob = c(1-p,p))
    
    return(a)
  }

  ggplotvor<-function(plotpoints,tit,wid = 5){
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
  
  #comienza el programa
  
  
  metropolisad_ben<-function(seed = 666, steps = 250, n = 100, Layers = 5,
                         RadiusA = 5/(2*pi), Ratio = 2.5, cyl_length = 20,
                         gamma_ad = 0.15, lambda_ad = 0.04, s0 = 1,
                         alpha = 1, beta = 100){
    
    
    #We define our variables
    
    RadiusB <- Ratio*RadiusA
    cyl_width_A <- 2*pi*RadiusA
    cyl_width_B <- 2*pi*RadiusB
    
    cyl_thickness <- RadiusB-RadiusA
    
    #We define the vertices of the apical plane
    
    xmin <- 0
    xmax <- cyl_width_A
    ymin <- 0
    ymax <- cyl_length
    
    r <- cyl_width_A/n #radius to make the moves
    Am <- ((RadiusA+RadiusB)*pi*cyl_length)/n
    
    rec <- list()
    rad <- list()
    
    for(k in 1:Layers){
      rad[[k]]<- RadiusA+(k-1)*(cyl_thickness/(Layers-1)) #the radius of the layer k
      rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
    }
    
    #Start, first iteration
    
    set.seed(seed)
    
    x1 <- runif(n,xmin,xmax)
    y1 <- runif(n,ymin,ymax)
    x <-c(x1,x1+cyl_width_A,x1+2*cyl_width_A)
    y <-c(y1,y1,y1)
    
    points <- vector(mode = "list", length = Layers)
    points <- lapply(1:Layers, function(i){data.frame(x = (rad[[i]]/rad[[1]])*x, y = y )})
    
    pointsinit <- points
    energytesel <- bending_tesellation_energy_N(points, Am, rec , rad,
                                        gamma_ad, lambda_ad, alpha, n, Layers, s0)
    
    #We create the variables to store the results
    
    energhist <- data.frame(iteration=numeric(steps), energy=numeric(steps))
    energhist[1,c(1,2)] <- c(0,energytesel)
    
    histpts <- vector(mode="list", length = Layers)
    
    for (i in 1:Layers) {
      histpts[[i]] <- data.frame(x = double(3*n*steps), y = double(3*n*steps), Frame = double(3*n*steps))
    }
    
    #Start of the loop
    
    for (j in 1:steps) {
      for(l in 1:n) {
        points2 <- bending_move_points(points, cyl_width_A, cyl_length, r,
                                       n, Layers, rad)
        
        energytesel2 <- bending_tesellation_energy_N(points2, Am, rec, rad,
                                           gamma_ad, lambda_ad, alpha, n,
                                           Lay = Layers, s0)
        c <- choice_metropolis(energytesel2-energytesel, beta)
        cond <- c==1
        if(cond){
          points <- points2
          energytesel <- energytesel2
        }
        gc()
      }
      for (i in 1:Layers) {
        histpts[[i]][(j*3*n+1):(j*3*n+3*n),c(1,2)] <- points[[i]]
        histpts[[i]][(j*3*n+1):(j*3*n+3*n),3] <- j+1
      }
      energhist[j+1,c(1,2)] <- c(j,energytesel)
      gc()
    }
    save(histpts, file = paste0("results_", i, ".Rds"))
    nu2 <- nu_sq(points = points, rec = rec, n = 100)
    return(list(histpts,energhist,nu2))
  }
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  results <- foreach(i=1000:1003,
                     .combine = rbind, .packages = "deldir") %dopar% {
    do.call(metropolisad_ben, list(seed = i, steps = 2, L=5, Ratio = 10))
  }
  
  stopCluster(cl)
  
  stopImplicitCluster()
  
  save(results, file = "results.Rds")