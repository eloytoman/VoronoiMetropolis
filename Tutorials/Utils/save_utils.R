library(deldir)
library(dplyr)


save_tessellation_layers <- function(pts,
                                     n = 100, RadiusA = 5/(2*pi), Ratio = 2.5,
                                     cyl_length = 20,
                                     Layers = 10, filename){
  RadiusB <- Ratio*RadiusA
  cyl_width_A <- 2*pi*RadiusA
  cyl_width_B <- 2*pi*RadiusB
  cyl_thickness <- RadiusB-RadiusA
  
  xmin <- 0
  xmax <- cyl_width_A
  ymin <- 0
  ymax <- cyl_length
  
  rec <- list()
  rad <- list()
  
  for(k in 1:Layers){
    rad[[k]]<- RadiusA+(k-1)*(cyl_thickness/(Layers-1)) #the radius of the layer k
    rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
  }
  
  df <- data.frame(id_cell= integer(), Layer = double(), radius = double(),
                   centroidx = double(), centroidy = double(),
                   n_vertices=integer(), vert1x = double(), vert2x = double(),
                   vert3x = double(), vert4x = double(), vert5x = double(),
                   vert6x = double(), vert7x = double(), vert8x = double(),
                   vert9x = double(), vert10x = double(), vert11x = double(),
                   vert1y = double(), vert2y = double(), vert3y = double(),
                   vert4y = double(), vert5y = double(), vert6y = double(), 
                   vert7y = double(), vert8y = double(), vert9y = double(), 
                   vert10y = double(), vert11y = double())
  
  for (j in 1:Layers) {
    tes <- deldir(pts$x*(rad[[j]]/rad[[1]]),pts$y,rw=rec[[j]])
    tiles <- tile.list(tes)[(n+1):(2*n)]
    lon<-length(df$id_cell)
    for (i in 1:n) {
      df[lon+i,1] <- tiles[[i]][[1]]
      df[lon+i,2] <- j
      df[lon+i,3] <- rad[[j]]
      df[lon+i,4] <- tiles[[i]][[2]][[1]]
      df[lon+i,5] <- tiles[[i]][[2]][[2]]
      df[lon+i,6] <- length(tiles[[i]][[3]])
      df[lon+i,7:(6+length(tiles[[i]][[3]]))] <- tiles[[i]][[3]]
      df[lon+i,18:(17+length(tiles[[i]][[4]]))] <- tiles[[i]][[4]]
    }
    gc()
    rm(tiles,tes)
  }
  
  write.table(df, file = filename, sep = ",", row.names = FALSE)
  return(df)
}

save_tessellation_layers_bending <- function(pts,
                                             n = 100, RadiusA = 5/(2*pi), Ratio = 2.5,
                                             cyl_length = 20,
                                             Layers=10, filename){
  RadiusB <- Ratio*RadiusA
  cyl_width_A <- 2*pi*RadiusA
  cyl_width_B <- 2*pi*RadiusB
  cyl_thickness <- RadiusB-RadiusA
  
  xmin <- 0
  xmax <- cyl_width_A
  ymin <- 0
  ymax <- cyl_length
  
  rec <- list()
  rad <- list()
  
  for(k in 1:Layers){
    rad[[k]]<- RadiusA+(k-1)*(cyl_thickness/(Layers-1)) #the radius of the layer k
    rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
  }
  
  df <- data.frame(id_cell = integer(), Layer = double(), Radius = double(),
                   Centroidx = double(), Centroidy = double(),
                   n_vertices = integer(), vert1x = double(), vert2x = double(),
                   vert3x = double(), vert4x = double(), vert5x = double(),
                   vert6x = double(), vert7x = double(), vert8x = double(),
                   vert9x = double(), vert10x = double(), vert11x = double(),
                   vert1y = double(), vert2y = double(), vert3y = double(),
                   vert4y = double(), vert5y = double(), vert6y = double(), 
                   vert7y = double(), vert8y = double(), vert9y = double(), 
                   vert10y = double(), vert11y = double())
  
  for (j in 1:Layers) {
    tes <- deldir(pts[[j]]$x, pts[[j]]$y, rw=rec[[j]])
    tiles <- tile.list(tes)[(n+1):(2*n)]
    lon <- length(df$Layer)
    for (i in 1:n) {
      df[lon+i,1] <- tiles[[i]]$ptNum-n
      df[lon+i,2] <- j
      df[lon+i,3] <- rad[[j]]
      df[lon+i,4] <- tiles[[i]][[2]][[1]]
      df[lon+i,5] <- tiles[[i]][[2]][[2]]
      df[lon+i,6] <- length(tiles[[i]][[3]])
      df[lon+i,7:(6+length(tiles[[i]][[3]]))] <- tiles[[i]][[3]]
      df[lon+i,18:(17+length(tiles[[i]][[4]]))] <- tiles[[i]][[4]]
    }
    gc()
    rm(tiles,tes)
  }
  
  # write.csv2(df, file = "output_databending.csv")
  write.table(df, file = filename, sep = ",", row.names = FALSE)
  return(df)
}
