library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)
library(tidyr)


energy_analysis <- function(points){
  
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Bending_energy = double(Lay),
                         Total_energy=double(Lay))
  
  for (i in 1:Lay) {
    tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
    tesener[i,c(2,3,4,6)] <- c(elener/Lay,
                               tenener/lay,
                               contener/Lay,
                               sum((areas-1)^2+(gam/2)*(perims^2)+
                                     lamad*perims)/(n*Lay))
    
  }
  
  tesener$Bending_energy[c(1,Lay)]<-0
  
  tesener$Bending_energy[2:(Lay-1)] <- sapply(2:(Lay-1), function(i){
    
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
    
    return(sum(alpha*((angles-pi)^2))/(Lay*n))
  })
  
  tesener$Total_energy <- tesener$Total_energy+tesener$Bending_energy
  return(tesener)
}

energy_analysis_nobend <- function(histpts, it=150, n=100, Lay=10){
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Total_energy=double(Lay))
  tesener2 <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                          Elastic_energy = double(Lay),
                          Tension_energy = double(Lay),
                          Contractile_energy = double(Lay),
                          Total_energy=double(Lay))
  
  points <- filter(histpts, Frame == 1)
  pointslast <- filter(histpts, Frame == it)
  
  for (i in 1:Lay) {
    tesel <- deldir(points$x*(rad[[i]]/rad[[1]]),points$y,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
    tesener[i,c(2,3,4,5)] <- c(elener/Lay,
                               tenener/Lay,
                               contener/Lay,
                               (elener+tenener+contener)/(Lay*n))
    
    tesellast <- deldir(pointslast$x*(rad[[i]]/rad[[1]]),pointslast$y,rw=rec[[i]])
    tilestlast <- tile.list(tesellast)[(n+1):(2*n)]
    
    perimsl <- (tilePerim(tilestlast)$perimeters)/sqrt(A0)
    areasl <- sapply(tilestlast,function(x){x$area/A0})
    
    elenerl <- (sum((areasl-1)^2)/n)
    tenenerl <- sum(lamad*perimsl)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contenerl <- sum((gam/2)*(perimsl^2))/n
    tesener2[i,c(2,3,4,5)] <- c(elenerl/Lay,
                                tenenerl/Lay,
                                contenerl/Lay,
                                (elenerl+tenenerl+contenerl)/(n*Lay))
  }
  
  p<-ggplot(tesener, aes(x = Layer))+
    geom_line(aes(y = Total_energy, colour = "Total energy"))+
    geom_line(aes(y=Tension_energy, colour = "Tension energy"))+
    geom_line(aes(y = Elastic_energy,colour = "Adhesion energy"))+
    geom_line(aes(y = Contractile_energy, colour = "Contractile energy"))+
    # geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
    labs(title = "Energy Analysis before executing the algorithm",
         x = "Layer ratio (R/Ra)", y = "Average energy of the cell",
         color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("Tension energy",
                                   "Adhesion energy",
                                   "Contractile energy",
                                   "Total energy"),
                        # "Bending energy"),
                        values = c("red",
                                   "blue",
                                   # "darkgreen",
                                   "orange",
                                   "purple"))
  
  show(p)
  
  p2<-ggplot(tesener2, aes(x = Layer))+
    geom_line(aes(y = Total_energy, colour = "Total energy"))+
    geom_line(aes(y=Tension_energy, colour = "Tension energy"))+
    geom_line(aes(y = Elastic_energy,colour = "Adhesion energy"))+
    geom_line(aes(y = Contractile_energy, colour = "Contractile energy"))+
    # geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
    labs(title = paste0("Energy Analysis after ",it, " iterations of the algorithm"),
         x = "Layer ratio (R/Ra)", y = "Average energy of the cell",
         color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("Tension energy",
                                   "Adhesion energy",
                                   "Contractile energy",
                                   "Total energy"),
                        # "Bending energy"),
                        values = c("red",
                                   "blue",
                                   # "darkgreen",
                                   "orange",
                                   "purple"))
  
  show(p2)
  
  
  return(list(tesener,tesener2))
}

energy_iteration <- function(pointsx, pointsy, n=100, Lay = 10){
  
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Bending_energy = double(Lay),
                         Total_energy=double(Lay))
  elasticener<-0
  tensionener<-0
  contractilener<-0
  totalener<-0
  
  for (i in 1:Lay) {
    tesel <- deldir(pointsx*(rad[[i]]/rad[[1]]),pointsy,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
    elasticener <- elasticener+elener
    tensionener <- tensionener+tenener
    contractilener <- contractilener+contener
    totalener<- totalener+elener+tenener+contener
  }
  elasticener <- elasticener/Lay
  tensionener <- tensionener/Lay
  contractilener <- contractilener/Lay
  totalener <- totalener/Lay
  
  return(c(elasticener,tensionener,contractilener,totalener))
}

energy_analisis_1sim <- function(histpts, it = 150, lay = 10, n =100){
  histener<-data.frame(it = 1:it, elen = double(it), tenen = double(it),
                       conten = double(it), toten = double(it))
  
  for (i in 1:it) {
    pts <- filter(histpts, Frame==i)
    histener[i,c(2,3,4,5)] <- energy_iteration(pts$x,pts$y, n = n, Lay = lay)
  }
  return(histener)
}

plot_energyss<-function(enerhist){
  
  p<-ggplot(enerhist, aes(x = it))+
    # geom_line(aes(y = toten, colour = "Total energy"))+
    geom_line(aes(y= tenen, colour = "Tension energy"))+
    geom_line(aes(y = elen, colour = "Adhesion energy"))+
    geom_line(aes(y = conten, colour = "Contractile energy"))+
    # geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
    labs(title = "Decomposition of system energies",
         x = "Iteration", y = "Average energy per cell",
         color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("Tension energy",
                                   "Adhesion energy",
                                   "Contractile energy"),
                        # "Total energy",
                        # "Bending energy"),
                        values = c("red",
                                   "blue",
                                   "darkgreen"))
  # "orange",
  # "purple"))
  
  show(p)
}

energy_layers_sim1 <- function(histpts,it = 150, Lay = 10, n = 100, rect = rec){
  
  tesener <- data.frame( iter = 1:it)
  
  for (i in 1:it) {
    pts<-filter(histpts,Frame==i)
    ener<-rep(0,6)
    k<-1
    layers_energy <- seq(1,Lay, by =2)
    for (j in layers_energy) {
      ener[[k]] <- energy_1_layer(pts$x*(rad[[j]]/rad[[1]]),pts$y,rect[[j]],
                                  i = j, Lay = Lay, n = n)
      k<-k+1
    }
    layers_len <- seq(2,length(layers_energy)+1,by = 1)
    
    tesener[i,layers_len] <- ener
  }
  colnames(tesener)<- c(c("iter"), paste0("layer_",layers_energy))
  return(tesener)
}

energy_1_layer <- function(pointsx, pointsy, rec, i, Lay =10, n=100){
  
  tesel <- deldir(pointsx,pointsy,rw=rec)
  tilest <- tile.list(tesel)[(n+1):(2*n)]
  
  perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
  areas <- sapply(tilest,function(x){x$area/A0})
  
  gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
  
  tesener<-sum((areas-1)^2+(gam/2)*(perims^2)+
                 lamad*perims)/(Lay*n)
  return(tesener)
}

plot_energy_decomp <- function(layers_hist) {
  df_largo <- pivot_longer(layers_hist, cols = starts_with("layer_"), names_to = "variable", values_to = "valor")
  
  ggplot(df_largo, aes(x = iter, y = valor, color = variable)) +
    geom_line() +
    labs(title="Energy relaxation by layer",x = "Iteration", y = "Energy", color = "Layer")
  
}
