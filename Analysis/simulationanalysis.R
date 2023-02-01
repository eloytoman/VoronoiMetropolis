library(deldir)
library(ggplot2)
library(ggvoronoi)
library(stats)
library(dplyr)
library(nls.multstart)


#This script takes the results after performing dinstinct simulations
#and makes the plots corresponding to the lewis law and the histogram of 
#sides of al the 400ths tessellations of every simulation

funaux<-function(p,rec = c(0,15,0,20), n = 100){
  
  #funaux returns a dataframe which first column is the edges of the cell and the
  # second column is the area of that cell
  
  tsl<-deldir(p$x,p$y,rw=rec)
  til<-tile.list(tsl)
  celledgearea<-data.frame(edges=integer(),area=double())
  for (i in 1:n) {
    celledgearea[i,c(1,2)]<-c(length(til[[i+n]]$x),til[[i+n]]$area/A0)
  }
  return (celledgearea)
}

funaux2sim<-function(ptsord, n = 100, ps = 100){
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  edgear<-data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a<-length(edgear$edges)
    edgear[(a+1):(a+n),c(1,2)]<-funaux(ptsord[((i-1)*(3*n)+1):(i*3*n),c(1,2)], n)
    edgear[(a+1):(a+n),3]<-i
  }
  return(edgear)
}

funaux2simDOUBLE<-function(ptsord, n = 100, ps = 100, Ratio = 2.5){
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  #This function returns a big dataframe, made up of one dataframe for each simulation
  # the df for one single simulation represents the edges and area of each of the cells.
  
  edgearA <- data.frame(edges=integer(), area=double(), Frame=integer())
  for (i in 1:ps) {
    a <- length(edgearA$edges)
    edgearA[(a+1):(a+n),c(1,2)] <- funaux(ptsord[((i-1)*(3*n)+1):(i*3*n),c(1,2)], rec=c(0,15,0,20), n)
    edgearA[(a+1):(a+n),3] <- i
  }
  
  ptsordB <- ptsord
  ptsordB$x <- ptsordB$x*Ratio
  
  edgearB <- data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a <- length(edgearB$edges)
    edgearB[(a+1):(a+n),c(1,2)] <- funaux(ptsordB[((i-1)*(3*n)+1):(i*3*n),c(1,2)], rec=c(0,15*Ratio,0,20))
    edgearB[(a+1):(a+n),3]<-i
  }
  return(list(edgearA,edgearB))
}

stationarylewis<-function(edgear){
  #First execute ord function (for results) and funaux2 with the data, then this function makes the plots.
  plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area))+
    geom_point() + xlab("Number of sides") + ylab("Relative area")+
    ggtitle("Relative area of the cells by sides for the system in equilibrium. Basal surface")+
    stat_summary(aes(y = area, group = 1), fun = mean, colour = "#00BFC4", geom = "line", group = 1)
  show(plotareaedges)
  
  mined<-min(edgear$edges)
  maxed<-max(edgear$edges)
  quant<-maxed-mined
  
  datahist<-data.frame(edges=numeric(),frec=numeric())
  
  for (i in min(edgear$Frame):max(edgear$Frame)){
    for (j in mined:maxed) {
      dat<-dplyr::filter(edgear, edges == j & Frame == i)
      pos<-(i-1)*quant+j-mined+1
      datahist[pos,c(1,2)]<-c(j, length(dat$edges))
    }
  }
  
  meandat <- data.frame(edges=numeric(), meanfrec=numeric())
  varcoefdat <- data.frame(edges=numeric(),varcoef=numeric())
  
  for (j in mined:maxed) {
    dat2<-dplyr::filter(datahist, edges==j)
    pos<-j-mined+1
    meandat[pos,c(1,2)]<-c(j,mean(dat2$frec))
    varcoefdat[pos,c(1,2)]<-c(j,var(dat2$frec)/mean(dat2$frec))
  }
  
  print(datahist)
  print(meandat)
  print(varcoefdat)
  
  histedges<-ggplot(meandat,aes(edges,meanfrec/sum(meanfrec)))+
    geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
    xlab("Number of edges")+ylab("Average frequency")+
    ggtitle("Average quantity of sides of the cells for system in equilibrium. Basal surface")+
    xlim(2,11)
  show(histedges)
}

ord <- function(results, iter = 150, sim=100){
  #With this function we extract the points for one iteration of every simulation
  #iter is the specific iteration of the algorithm that we extract to make the analysis
  
  ptsord<-data.frame(x=numeric(300*sim),y=numeric(300*sim))
  for (i in 0:(sim-1)) {
    ptsord[(i*300+1):((i+1)*300),c(1,2)]<-results[[i+1]][results[[i+1]]$Frame==iter,c(1,2)]
  }
  return(ptsord)
}

regnls<-function(energh){
  
  ##This functions makes a single regression for one iteration with the function
  #nls
  
  x<-unlist(lapply(energh$iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls(y~I(a+b*(1-exp(-x/c))),start = list(a=2,b=0.3,c=0.1))
  #plot(x,y)
  #lines(x,predict(m),col="red",lwd=3)
  summary(m)
}

regnls2<-function(energh, n=100){
  
  #This functions makes a single regression for one iteration with the function
  #nls_multistart (which is much better for all cases)
  
  x<-unlist(lapply(energh$iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls_multstart(y~I(a+b*(1-exp(-x/c))),
                   iter = 500,
                   start_lower = list(a=0,b=-5,c=5),
                   start_upper = list(a=5,b=5,c=30))
  #plot(x,y)
  #lines(x,predict(m),col="red",lwd=3)
  summary(m)
  return(summary(m)[["coefficients"]][c(1,2,3)])
}

adjsim<-function(results,nsim=100,it=150){
  
  #This function makes a regression of the energy of the results for each iteration,
  # and then it makes the average.
  
  coefest<-data.frame(a=double(nsim),b=double(nsim),c=double(nsim))
  for (i in 1:nsim) {
    coefest[i,c(1,2,3)]<-regnls2(results[[i+nsim]])
  }
  a<-mean(coefest$a)
  b<-mean(coefest$b)
  c<-mean(coefest$c)
  adj<-function(x) a+b*(1-exp(-x/c))
  adj_data<-data.frame(x=1:it,y=adj(1:it))
  ploten<-ggplot(data = adj_data, aes(x=x,y=y))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Average energy relaxation of the system. 100 simulations")
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

scutoids_percr<-function(histpts, rect1, rect2, n = 100, it=150){
  #Computes the evolution of the percentage of escutoids in every iteration of the algorithm
  perc<-data.frame(it=integer(it), perc_sc=double(it))
  for (i in 1:it) {
    pointsAx<-filter(histpts, Frame == i)$x
    pointsAy<-filter(histpts, Frame == i)$y
    pointsBx <- pointsAx*2.5
    pointsBy <- pointsAy*2.5
    tslA<-deldir(pointsAx,pointsAy,rw=rect1)
    tilA<-tile.list(tslA)[(n+1):(2*n)]
    tslB<-deldir(pointsBx,pointsBy,rw=rect2)
    tilB<-tile.list(tslB)[(n+1):(2*n)]
    cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
    for (j in 1:length(tilA)) {
      cellsdf[j,c(1,2)]<-c(length(tilA[[j]]$x),length(tilB[[j]]$x))
    }
    percen<-length(filter(cellsdf, edgesA==edgesB)[[1]])
    perc[i,c(1,2)]<-c(i,percen)
  }
  # ploten<-ggplot(perc,aes(x=it,y=perc_sc))+
  #   geom_line(colour="#F8766D")+
  #   xlab("Iteration of the algorithm")+
  #   ylab("Percentage of escutoids")+
  #   ggtitle("Evolution of the percentage of escutoids of the system")
  # show(ploten)
  return(perc)
}

scutoids_percr_simulations<-function(results, rect1, rect2, n=100, sim=100, it=150){
  perc<-data.frame(it=1:it, percen=rep(0,it))
  for (j in 1:sim) {
    percen<-scutoids_percr(results[[j]], rect1, rect2, n, it)$perc_sc
    perc$percen<-perc$percen + percen
  }
  perc$percen<-perc$percen/sim
  ploten<-ggplot(perc,aes(x=it,y=percen))+
      geom_line(colour="#F8766D")+
      xlab("Iteration of the algorithm")+
      ylab("Percentage of escutoids")+
      ggtitle("Evolution of the average percentage of escutoids. 100 simulations")
    show(ploten)
  
  return(perc)
}

scutoids_percr_simulations(results,rec1,rec2)

scutoids_prep <- function(pointsAx,pointsAy,pointsBx,pointsBy,rect1,rect2,
                          n = 100){
  
  #This function returns a dataframe which counts the number of edges in the 
  # apical and basal surface of each cell, and counts the cells that has the same
  # number of edges in both surfaces.
  
  tslA <- deldir(pointsAx,pointsAy,rw=rect1)
  tilA <- tile.list(tslA)[(n+1):(2*n)]
  tslB <- deldir(pointsBx,pointsBy,rw=rect2)
  tilB <- tile.list(tslB)[(n+1):(2*n)]
  cellsdf <- data.frame(edgesA=integer(),edgesB=integer())
  for (i in 1:n) {
    cellsdf[i,c(1,2)] <- c(length(tilA[[i]]$x),length(tilB[[i]]$x))
  }
  countdf <- cellsdf %>%
    count(edgesA,edgesB)
  return(countdf)
}

scutoids_analysis_stationary <- function(histpts, rect1, rect2, n = 100){
  lon <- 50 #how many iterations we want to have
  histdf_count <- data.frame(edgesA=integer(lon*n),edgesB=integer(lon*n),count=integer(lon*n))
  for (i in 1:lon) {
    histdf_count[((i-1)*n+1):(i*n)]<-
      scutoids_prep(dplyr::filter(histpts, Frame==250+i-1)$x,
                    dplyr::filter(histpts, Frame==250+i-1)$y,
                    (Radius2/Radius)*dplyr::filter(histpts, Frame==250+i-1)$x,
                    dplyr::filter(histpts, Frame==250+i-1)$y,
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

scutoids_analysis_simulations <- function(results, Ratio = 2.5, rect1 = rec1, rect2 = rec2,
                                          n = 100, sim = 100, it = 150){
  
  #sim is how many simulations we have
  #it is the iteration that we have
  
  
  histdf_count <- data.frame(edgesA=integer(),
                             edgesB=integer(),
                             count=integer())

  for (i in 1:sim) {
    #each simulation have a different distribution of sides, so we have to adapt
    
    a <- length(histdf_count[[1]])
    ptsx <- dplyr::filter(results[[i]], Frame==it)$x
    ptsy <- dplyr::filter(results[[i]], Frame==it)$y
    df<-scutoids_prep(ptsx, ptsy,
                      Ratio*ptsx, ptsy,
                      rect1,rect2)
    len<-length(df[[1]])
    histdf_count[(a+1):(a+len),c(1,2,3)]<-df
  }
  
  #Now, to compute the average, we have to take into account the 1000 simulations
  
  histdf_avcount <- data.frame(edgesA = double(),
                               edgesB = double(),
                               avg_count = double())
  for (edA in min(histdf_count$edgesA):max(histdf_count$edgesA)) {
    for (edB in min(histdf_count$edgesB):max(histdf_count$edgesB)) {
      a <- length(histdf_avcount$edgesA)
      histdf_avcount[a+1,c(1,2,3)] <-
        c(edA, edB, (sum(dplyr::filter(histdf_count, edgesA == edA & edgesB == edB)$count)/100))
    }
  }
  histdf_avcount <- dplyr::filter(histdf_avcount, avg_count >= 0.1)
  
  scutoidsplot<-ggplot(histdf_avcount, aes(x = edgesA, y = edgesB, label=avg_count))+
    geom_count(shape = "square", aes(color= avg_count))+
    xlab("Average of edges on apical surface")+ylab("Average of edges on basal surface")+
    ggtitle("Average polygon class of apical and basal surfaces")+
    guides(colour = "colorbar", size = "none")+
    scale_size_area(max_size = 30)+
    geom_label()+
    scale_fill_gradient(low = "light blue", high = "deepskyblue")+
    xlim(3.5,8.5)+
    ylim(3.5,8.5)
    # xlim(min(histdf_avcount$edgesA)-0.5, max(histdf_avcount$edgesA)-0.5)+
    # ylim(min(histdf_avcount$edgesB)-0.5, max(histdf_avcount$edgesB)-0.5)
  show(scutoidsplot)
  return(histdf_avcount)
}

resord<-ord(results)
edgearsim<-funaux2simDOUBLE(resord)
stationarylewis(edgearsim[[1]][1:10000,c(1,2,3)])
stationarylewis(edgearsim[[2]][1:10000,c(1,2,3)])
coef<-adjsim(results)

xmin<-0
xmax<-5
ymin<-0
ymax<-20
rec1 <- c(xmin,xmin+3*xmax,ymin,ymax)
rec2 <- c(xmin,xmin+3*(xmax*2.5),ymin,ymax)



scutoids_prep(points$x,points$y,(2.5)*points$x,points$y,rect1= rec1,rect2=rec2)
scutoids_analysis_oneiter(points$x,points$y,(2.5)*points$x,points$y,rect1= rec1,rect2=rec2)
df_scutoid <- scutoids_analysis_simulations(results)

save_tessellation <- function(pts, rec = c(xmin,xmin+3*xmax,ymin,ymax),
                              n = 100, radius = 5/(2*pi)){
  b <- 1
  tes <- deldir(pts$x,pts$y,rw=rec)
  tiles <- tile.list(tes)[(n+1):(2*n)]
  df <- data.frame(cell_id = double(), radius = double(),
                   centroidx = double(), centroidy = double(),
                   n_vertices=integer(), vert1x = double(), vert2x = double(),
                   vert3x = double(), vert4x = double(), vert5x = double(),
                   vert6x = double(), vert7x = double(), vert8x = double(),
                   vert9x = double(), vert10x = double(), vert11x = double(),
                   vert1y = double(), vert2y = double(), vert3y = double(),
                   vert4y = double(), vert5y = double(), vert6y = double(), 
                   vert7y = double(), vert8y = double(), vert9y = double(), 
                   vert10y = double(), vert11y = double())
  for (i in 1:n) {
    df[i,1] <- tiles[[i]][[1]]
    df[i,2] <- radius
    df[i,3] <- tiles[[i]][[2]][[1]]
    df[i,4] <- tiles[[i]][[2]][[2]]
    df[i,5] <- length(tiles[[i]][[3]])
    df[i,6:(5+length(tiles[[i]][[3]]))] <- tiles[[i]][[3]]
    df[i,17:(16+length(tiles[[i]][[4]]))] <- tiles[[i]][[4]]
  }
  write.csv2(df, file = "output_data2.csv")
  write.csv(df,file = "output_data.csv")
  return(df)
}

# metropolisad_ben<-function(seed = 666, steps = 250, n = 100, Layers = 5,
#                            RadiusA = 5/(2*pi), Ratio = 2.5, cyl_length = 20,
#                            gamma_ad = 0.15, lambda_ad = 0.04, s0 = 1,
#                            alpha = 1, beta = 100,
# ){




save_tessellation_Layers <- function(pts,
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




save_tessellation_Layers(points, filename = "output_data.csv")





save_tessellation_Layers_Bending <- function(pts,
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



dframe=save_tessellation_Layers_Bending(pts = points, Layers = 12,
                                        filename = "output_databending_alpha2_150it.csv")




x1 <- runif(n,xmin,xmax)
y1 <- runif(n,ymin,ymax)

x <- c(x1,x1+cyl_width,x1+2*cyl_width)
y <- c(y1,y1,y1)

#We create a list of dataframes which will store the points in each layer.
pointsinit <- vector(mode = "list", length = 12)
pointsinit <- lapply(1:L, function(i){data.frame(x=(rad[[i]]/Radius)*x,y=y)})


pointsinit <- filter(results[[1]], Frame ==1)
pointsfin <- filter(results[[1]], Frame==550)
save_tessellation_Layers(pointsinit, filename = "nobend_f1.csv")
save_tessellation_Layers(pointsfin, filename = "nobend_f150.csv")

alphas<-c(0.05,0.2,0.5,1,2,5,10,50)

for (i in 1:8) {
  for (j in 1:Lay) {
    points[[j]]<- filter(results[[i]][[j]], Frame==300)
  }
  tit<-paste0("data_alpha_",alphas[[i]],".csv")
  save_tessellation_Layers_Bending(points, filename = tit)
}

for (j in 1:Lay) {
  points[[j]]<- filter(results[[6]][[j]], Frame==550)
}

save_tessellation_Layers_Bending(points, filename = "cells05.csv")


scutoids_perc_cells_1layer <- function(points1x,points1y,points2x,points2y, rec1,rec2, n = 100){
  
  tslA<-deldir(points1x,points1y,rw=rec1)
  tilA<-tile.list(tslA)[(n+1):(2*n)]
  tslB<-deldir(points2x,points2y,rw=rec2)
  tilB<-tile.list(tslB)[(n+1):(2*n)]
  
  cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
  countdf<-data.frame(cell=1:100, intercalations=rep(0,100))
  
  for (i in 1:length(tilA)) {
    cellsdf[i,c(1,2)]<-c(length(tilA[[i]]$x),length(tilB[[i]]$x))
  }
  for (i in 1:100) {
    if(cellsdf$edgesA[[i]]!=cellsdf$edgesB[[i]]){countdf$intercalations[[i]]=1}
  }
  return(countdf)
}


scutoids_perc_cells_sim_nobend <- function(points, rec, rad, Lay=10){
  
  countscut<-data.frame(cell=1:100, intercalations=rep(0,100))
  for(j in 1:(Lay-1)){
    countscut$intercalations <- countscut$intercalation +
      scutoids_perc_cells_1layer(points$x*(rad[[j]]/rad[[1]]),points$y,
                                 points$x*(rad[[j+1]]/rad[[1]]),points$y,
                                 rec[[j]],rec[[j+1]])
  }
  return(countscut)
}

points<-filter(histpts, Frame == 300)
sc_pc <- scutoids_perc_cells_sim_nobend(points, rec, rad)
perc_sc <-length(filter(sc_pc), intercalations>0)
print(perc_sc)


