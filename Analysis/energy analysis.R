library(deldir)
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)


# This functions computes the energy of a certain tessellation, given the
# points that generate it. It makes a distinction between the different types 
# of energies, (bending, elastic, tension and contractile). So this functions 
# makes the energy decomposition in a given instant of the algorithm for every 
# layer of the cylinder. The output is a dataframe which has the info of every
# type of energy for all layers in the given instant.

energy_analysis <- function(points, Lay){
  
  
  #First we initialize the dataframe that store the results.
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Bending_energy = double(Lay),
                         Total_energy=double(Lay))
  
  
  # Now we apply the same computation that we do when computing the energy 
  # during the algorithm, but storing every energy separately. To see more 
  # details, see the function "tesellation_energy_N" in N_cylinder_algorithm.R 
  # script
  
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
  
  #Same with the bending energy
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


# Example of use of the function above for the results of the parallel algorithm
# with bending, specifically, we use the 4th simulation and we take the
# iteration number 550 (when the algorithm has reached equilibrium). We 
# initialize the necessary constants and then we load the data from results.

gamad <- 0.15
lamad <- 0.04
alpha<-0.05

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
Lay <- 10  #Layers
lay <- Lay
r <- cyl_width/n 
bet <- 10
steps <- 550

A0 <- (((Radius2+Radius)/2)*2*pi*cyl_length)/n

rec <- list()
rad <- list()
for(k in 1:Lay){
  rad[[k]]<- Radius+(k-1)*(cyl_thickness/(Lay-1)) #the radius of the layer k
  rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
}

points <- list(Lay)
for (j in 1:Lay) {
  points[[j]]<-filter(results[[4]][[j]],Frame==550)
}

ener<-energy_analysis(points)

p<-ggplot(ener, aes(x = Layer))+
  geom_line(aes(y = Total_energy, colour = "Total energy"))+
  geom_line(aes(y = Tension_energy, colour = "Tension energy"))+
  geom_line(aes(y = Elastic_energy,colour = "Adhesion energy"))+
  geom_line(aes(y = Contractile_energy, colour = "Contractile energy"))+
  geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
  labs(title = "Energy Analysis after 550 iterations of the algorithm with Bending. Omega=0.05",
       x = "Layer Ratio (s=R/Ra)", y = "Average energy per cell",
       color = "Energy type") +
  scale_colour_manual("",
                      breaks = c("Tension energy",
                                 "Adhesion energy",
                                 "Contractile energy",
                                 "Bending energy",
                                 "Total energy"),
                      values = c("red",
                                 "blue",
                                 "darkgreen",
                                 "orange",
                                 "purple"))


show(p)






# Next function makes the energy decomposition in the first iteration and in 
# the last iteration of the algorithm, so we can compare the initial random
# distribution with the energy decomposition in an equilibrium state.
# In this case, we
# The function takes as input the histpts dataframe, which is the result of the
# N_cylinder_algorithm. This version does not take into account the bending 
# energy. The output are two graphics that show the energy decomposition 
# by layers at the beginning and at the end of the algorithm.

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
    labs(title = "Energy Analysis after 150 iterations of the algorithm",
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



# This function computes the value of every energy type for 1 iteration. The 
# objective is measure in a single iteration the values of every energy
# (tension, contractile, elastic and total energy).
# This one does not support bending energy.
# The difference with the first one is that the first one computes every energy
# for every layer in one iteration, and this one computes every energy for all
# layers together in one iteration (we sum energies of all layers)

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



# Given one whole simulation (we perform the algorithm twice and we store the 
# results in the dataframe of points histpts), this function stores the energy
# decomposition of every iteration (with the function above "energy_iteration").
# Then it returns the dataframe with the decomposition of energies for every
# iteration, having the historic of energy relaxation for every specific energy.
energy_analisis_1sim <- function(histpts, it = 150){
  histener<-data.frame(it = 1:it, elen = double(it), tenen = double(it),
                       conten = double(it), toten = double(it))
  
  for (i in 1:it) {
    pts <- filter(histpts, Frame==i)
    histener[i,c(2,3,4,5)] <- energy_iteration(pts$x,pts$y)
  }
  return(histener)
}

# This functions takes as input the output of the function above 
# "energy_analisis_1sim" and plot the energy relaxation for every specific
# energy
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

# Example of use of the functions above with hist1, the output of the script
# N_cylinder_algorithm.R
enerhist <- energy_analisis_1sim(hist1)
plot_energyss(enerhist)


#This function computes the total energy of every layer in one iteration
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

# This function takes as input the historic of points (output of 
# N_cylinder_algorithm) and returns the energy relaxation for every layer along
# the execution of the algorithm. In this case we chose only 6 layers, because
# the energy in intermediate layers is similar.
energy_layers_sim1 <- function(histpts,it = 150){
  
    tesener <- data.frame( iter = 1:it,
                           enlay1 = double(it),
                           enlay3 = double(it),
                           enlay5 = double(it),
                           enlay6 = double(it),
                           enlay8 = double(it),
                           enlay10 = double(it))
    
    for (i in 1:it) {
      pts<-filter(histpts,Frame==i)
      ener<-rep(0,6)
      k<-1
      for (j in c(1,3,5,6,8,10)) {
        ener[[k]] <- energy_1_layer(pts$x*(rad[[j]]/rad[[1]]),pts$y,rec[[j]], i = j)
        k<-k+1
      }
      tesener[i,c(2,3,4,5,6,7)] <- ener
    }
    return(tesener)
}

# This function plots the results of the function above, so we can visualize 
# the energy relaxation by layers
plot_energy_decomp<-function(tesenerdec){
  
ploten <- ggplot(tesenerdec, aes(x = iter))+
    geom_line(aes(y = enlay1, colour = "Apical layer (first)"))+
    geom_line(aes(y = enlay3, colour = "Third layer"))+
    geom_line(aes(y = enlay5, colour = "Fifth layer"))+
    geom_line(aes(y = enlay6, colour = "Sixth layer"))+
    geom_line(aes(y = enlay8, colour = "Eight layer"))+
    geom_line(aes(y = enlay10, colour = "Basal layer (tenth)"))+
    labs(title = "Decomposition of system energies by layer",
       x = "Iteration", y = "Average energy per cell",
       color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("Apical layer (first)",
                                   "Third layer",
                                   "Fifth layer",
                                   "Sixth layer",
                                   "Eight layer",
                                   "Basal layer (tenth)"),
                        values = c("red",
                                   "darkcyan",
                                   "magenta",
                                   "orange",
                                   "green",
                                   "blue"))

  show(ploten)
}





# The next function is used after performing a search of the best parameters
# as input for the algorithm, particulary changing the values of Omega.
# The function takes as input the results of the parallel algorithm 
# without bending (performed with different values of the Omega parameter) 
# and plots the energy relaxations for every simulation (every value of the 
# parameter)
plotenergies<-function(results,it=551){
  en<-data.frame(iteration=1:it, en1=results[[9]]$energy/100,
                 en2=results[[10]]$energy/100, en3=results[[11]]$energy/100,
                 en4=results[[12]]$energy/100, en5=results[[13]]$energy/100,
                 en6=results[[14]]$energy/100, en7=results[[15]]$energy/100,
                 en8=results[[16]]$energy/100)
  p<-ggplot(en, aes(x = iteration))+
    # geom_line(aes(y = Total_energy, colour = "Total energy"))+
    geom_line(aes(y = en1, colour = "0.001"))+
    geom_line(aes(y = en2, colour = "0.005"))+
    geom_line(aes(y = en3, colour = "0.01"))+
    geom_line(aes(y = en4, colour = "0.05"))+
    geom_line(aes(y = en5, colour = "0.1"))+
    geom_line(aes(y = en6, colour = "0.5"))+
    geom_line(aes(y = en7, colour = "1"))+
    geom_line(aes(y = en8, colour = "2"))+
    labs(title = "Energy Relaxation. Algorithm with Bending",
         x = "Iteration", y = "Average energy of the cells",
         color = "Value of Omega") +
    scale_colour_manual("",
                        breaks = c("0.001",
                                   "0.005",
                                   "0.01",
                                   "0.05",
                                   "0.1",
                                   "0.5",
                                   "1",
                                   "2"),
                        values = c("red",
                                   "blue",
                                   "magenta",
                                   "orange",
                                   "cyan",
                                   "darkgreen",
                                   "purple",
                                   "black"))
  p <- p + guides(colour=guide_legend(title="Values of Omega"))
  
  show(p)
  
}
plotenergies(results = results)



### BENDING FUNCTIONS



# This function computes the total energy of one single layer for the algorithm
# with bending. Here, "histpts" is the ouput of the bending algorithm,
# "rect" is the border of the cylinder in the layer we want to compute the energy,
# "rad" is the radius of the layer, "i" is the layer and "it" is the iteration of the 
# algorithm
energy_1_layer_bending <- function(histpts, rect, rad, i, it, Lay =10, n=100){
  
  points<-list(10)
  
  for (j in 1:Lay) {
    points[[j]]<- filter(histpts[[j]], Frame==it)
  }
  
  tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rect)
  tilest <- tile.list(tesel)[(n+1):(2*n)]
  
  perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
  areas <- sapply(tilest,function(x){x$area/A0})
  gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
  
  tesener<-sum((areas-1)^2+(gam/2)*(perims^2)+
                 lamad*perims)/(Lay*n)
  
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
  bendener <- sum(alpha*((angles-pi)^2))/(Lay*n)
  
  return(tesener+bendener)
}

# This function stores decomposes the energy relaxation by layers,
# so we have specifically the energy relaxation in every layer. As input, it 
# takes the histpts list of dataframes which is output of the algorithm with
# bending
bending_energy_layers_sim1 <- function(histpts,it =550){
  
  tesener <- data.frame( iter = 1:it,
                         enlay1 = double(it),
                         enlay2 = double(it),
                         enlay3 = double(it),
                         enlay4 = double(it),
                         enlay5 = double(it),
                         enlay6 = double(it),
                         enlay7 = double(it),
                         enlay8 = double(it),
                         enlay9 = double(it),
                         enlay10 = double(it))
  
  for (i in 2:it) {
    ener<-rep(0,6)
    k<-1
    for (j in c(2,3,4,5,6,7,8,9)) {
      ener[[k]] <- energy_1_layer_bending(histpts, rec[[j]], rad, i = j, it = i)
      k<-k+1
    }
    tesener[i,c(3,4,5,6,7,8,9,10)] <- ener
    points1 <- filter(histpts[[1]],Frame==i)
    points10 <- filter(histpts[[10]],Frame==i)
    tesener[i,2]<-energy_1_layer(points1$x,points1$y, rec[[1]],1)
    tesener[i,11]<-energy_1_layer(points10$x,points10$y,rec[[10]],10)
  }
  return(tesener)
}


# This function plots the energy relaxations for every layer. It takes as input
# the output of the previous function.
plot_energy_decomp_bending<-function(tesenerdec){
  
  ploten <- ggplot(tesenerdec, aes(x = iter))+
    geom_line(aes(y = enlay1, colour = "First layer"))+
    geom_line(aes(y = enlay2, colour = "Second layer"))+
    geom_line(aes(y = enlay3, colour = "Third layer"))+
    geom_line(aes(y = enlay4, colour = "Fourth layer"))+
    geom_line(aes(y = enlay5, colour = "Fifth layer"))+
    geom_line(aes(y = enlay6, colour = "Sixth layer"))+
    geom_line(aes(y = enlay7, colour = "Seventh layer"))+
    geom_line(aes(y = enlay8, colour = "Eight layer"))+
    geom_line(aes(y = enlay9, colour = "Ninth layer"))+
    geom_line(aes(y = enlay10, colour = "Tenth layer"))+
    labs(title = "Evolution of system energies by layer",
         x = "Iteration", y = "Average energy per cell",
         color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("First layer",
                                   "Second layer",
                                   "Third layer",
                                   "Fourth layer",
                                   "Fifth layer",
                                   "Sixth layer",
                                   "Seventh layer",
                                   "Eight layer",
                                   "Ninth layer",
                                   "Tenth layer"),
                        values = c("cyan",
                                   "darkgreen",
                                   "red",
                                   "darkcyan",
                                   "pink",
                                   "magenta",
                                   "orange",
                                   "green",
                                   "blue",
                                   "black"))
  
  show(ploten)
}


#Example of using the previous functions with one simulation performed with 
# the parallel algorithm (the fourth simulation)

hist1<-results[[4]]
tes<-bending_energy_layers_sim1(hist1)

plot_energy_decomp_bending(tes)



# This computes the energy decomposition by energy type in a single iteration
# and returns the energy of every type in the iteration. It does not make 
# distinctions by layer (all are sumed up)
energy_iteration_bending <- function(histpts, it, n=100, Lay = 10){
  
  points<-list(Lay)
  for (j in 1:Lay) {
    points[[j]]<- filter(histpts[[j]], Frame==it)
  }
  
  elasticener<-0
  tensionener<-0
  contractilener<-0
  totalener<-0
  bendener<-0
  
  for (i in 1:Lay) {
    tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
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


  for (i in 2:(Lay-1)) {

    angles <- numeric(n)

    for (j in 1:n) {
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

      v<-pmin(pmax(((vec1%*%vec2)[1,1])/
                     (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
      angles[j] <- acos(v)
    }
   bendener <- bendener+ sum(alpha*((angles-pi)^2))/(Lay*n)
  }
  totalener <- totalener+bendener
  
  return(c(elasticener,tensionener,contractilener,bendener,totalener))
}

# This function computes the energy relaxation for every energy type
energy_analisis_1sim_bending <- function(histpts, it = 300){
  histener<-data.frame(it = 1:it, elen = double(it), tenen = double(it),
                       conten = double(it), benden = double(it),
                       toten = double(it))
  
  for (i in 3:it) {
    histener[i,c(2,3,4,5,6)] <- energy_iteration_bending(histpts, it = i)
  }
  return(histener)
}

# This function plots the results from the previous function
plot_energyss_bending<-function(enerhist){
  
  p<-ggplot(enerhist, aes(x = it))+
    # geom_line(aes(y = toten, colour = "Total energy"))+
    geom_line(aes(y= tenen, colour = "Tension energy"))+
    geom_line(aes(y = elen, colour = "Adhesion energy"))+
    geom_line(aes(y = conten, colour = "Contractile energy"))+
    geom_line(aes(y = benden, colour = "Bending energy"))+
    labs(title = "Decomposition of system energies",
         x = "Iteration", y = "Average energy per cell",
         color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("Tension energy",
                                   "Adhesion energy",
                                   "Contractile energy",
                                   # "Total energy",
                                   "Bending energy"),
                        values = c("red",
                                   "blue",
                                   "darkgreen",
                                   "orange"))
  # "orange",
  # "purple"))
  
  show(p)
}


# Example with the data that we have from before

enerben <- energy_analisis_1sim_bending(hist1, 550)
plot_energyss_bending(enerben)



# Returns the total energy with bending of a certain tessellation (this is the 
# function used in the algorithm). The input is a list of dataframes, 
# each dataframe being the points in the iteration for a certain layer (each
# layer is stored in each element of the list)
bending_tesellation_energy_N <- function(points){
  
  tesener<-numeric(L)
  
  tesener<-sapply(1:L,function(i) {
    tesel<-deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perims<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas<-sapply(tilest,function(x){x$area/A0})
    gam<-gamma_ad*exp((1-(rad[[i]]/rad[[1]]))/s0)
    sum((areas-1)^2+(gam/2)*(perims^2)+lambda_ad*perims)
  })
  
  bendener <- sapply(2:(L-1), function(i){
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
      v <- pmin(pmax(((vec1%*%vec2)[1,1])/
                       (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
      ang <- acos(v)
      return(ang)
    })
    return(sum(alpha*((angles-pi)^2)))
  })
  
  return((sum(tesener)+sum(bendener))/L)
}


# Example of using the the function above to compute the energy in every 
# iteration after performing the algorithm with bending (recall that this
# algorithm have as output a list of dataframes, every dataframe being the
# energy relaxation for a layer)
enertes <- data.frame(iteration=3:550, energy= double(548))
for(k in 3:550){
  points<-list(10)
  for (j in 1:10) {
    points[[j]]<-filter(hist1[[j]],Frame==k)
  }
  enertes[k,2]<-bending_tesellation_energy_N(points)
}
plot_energy(enertes)
gamma_ad<-gamad
lambda_ad<-lamad
s0<-1
