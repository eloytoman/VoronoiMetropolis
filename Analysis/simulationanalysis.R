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

funaux2sim<-function(ptsord, n = 100, ps = 1000){
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

funaux2simDOUBLE<-function(ptsord, n = 100, ps = 1000, Ratio = 2.5){
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  #This function returns a big dataframe, made up of one dataframe for each simulation
  # the df for one single simulation represents the edges and area of each of the cells.
  
  edgearA<-data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a<-length(edgearA$edges)
    edgearA[(a+1):(a+n),c(1,2)]<-funaux(ptsord[((i-1)*(3*n)+1):(i*3*n),c(1,2)], rec=c(0,15,0,20), n)
    edgearA[(a+1):(a+n),3]<-i
  }
  
  ptsordB<-ptsord
  ptsordB$x<-ptsordB$x*Ratio
  
  edgearB<-data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a<-length(edgearB$edges)
    edgearB[(a+1):(a+n),c(1,2)]<-funaux(ptsordB[((i-1)*(3*n)+1):(i*3*n),c(1,2)], rec=c(0,15*Ratio,0,20))
    edgearB[(a+1):(a+n),3]<-i
  }
  
  return(list(edgearA,edgearB))
}

stationarylewis<-function(edgear){
  #First execute ord function (for results) and funaux2 with the data, then this function makes the plots.
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
      dat<-dplyr::filter(edgear, edges == j & Frame == i)
      pos<-(i-1)*quant+j-mined+1
      datahist[pos,c(1,2)]<-c(j, length(dat$edges))
    }
  }
  meandat<-data.frame(edges=numeric(), meanfrec=numeric())
  
  for (j in mined:maxed) {
    dat2<-dplyr::filter(datahist, edges==j)
    pos<-j-mined+1
    meandat[pos,c(1,2)]<-c(j,mean(dat2$frec))
  }
  print(datahist)
  print(meandat)
  
  histedges<-ggplot(meandat,aes(edges,meanfrec/sum(meanfrec)))+
    geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
    xlab("Number of edges")+ylab("Average frequency")+
    ggtitle("Average quantity of sides of the cells for the tessellations in stationary state")+
    xlim(2,11)
  show(histedges)
}

ord<-function(results, iter = 200){
  #With this function we extract the points for one iteration of every simulation
  #iter is the specific iteration of the algorithm that we extract to make the analysis
  
  ptsord<-data.frame(x=numeric(300*1000),y=numeric(300*1000))
  for (i in 0:999) {
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
  y<-unlist(lapply(energh$energy,as.numeric))/n
  m<-nls_multstart(y~I(a+b*(1-exp(-x/c))),
                   iter = 500,
                   start_lower = list(a=0,b=-5,c=5),
                   start_upper = list(a=5,b=5,c=30))
  #plot(x,y)
  #lines(x,predict(m),col="red",lwd=3)
  summary(m)
  return(summary(m)[["coefficients"]][c(1,2,3)])
}

adjsim<-function(results){
  
  #This function makes a regression of the energy of the results for each iteration,
  # and then it makes the average.
  
  coefest<-data.frame(a=double(1000),b=double(1000),c=double(1000))
  for (i in 1:1000) {
    coefest[i,c(1,2,3)]<-regnls2(results[[i+1001]])
  }
  a<-mean(coefest$a)
  b<-mean(coefest$b)
  c<-mean(coefest$c)
  adj<-function(x) a+b*(1-exp(-x/c))
  adj_data<-data.frame(x=1:200,y=adj(1:200))
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
  for (i in 1:length(tilA)) {
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


scutoids_analysis_simulations <- function(results, Ratio = 2.5, rect1 = rec, rect2 = rec2,
                                          n = 100, sim = 1000){
  
  #sim is how many simulations we have
  #it is the iteration that we have
  
  it <- 200
  
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
        c(edA, edB, (sum(dplyr::filter(histdf_count, edgesA == edA & edgesB == edB)$count)/1000))
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
    xlim(min(histdf_avcount$edgesA)-0.5, max(histdf_avcount$edgesA)-0.5)+
    ylim(min(histdf_avcount$edgesB)-0.5, max(histdf_avcount$edgesB)-0.5)
  show(scutoidsplot)
  return(histdf_avcount)
}


resord<-ord(results)
edgearsim<-funaux2simDOUBLE(resord)
stationarylewis(edgearsim[[1]][1:1000,c(1,2,3)])
stationarylewis(edgearsim[[2]][1:1000,c(1,2,3)])
coef<-adjsim(results)

scutoids_analysis_oneiter(points$x,points$y,(Radius2/Radius)*points$x,points$y,rec,rec2)
df_scutoid <- scutoids_analysis_simulations(results)

save_tessellation <- function(pts, rec = c(xmin,xmin+3*cyl_width,ymin,ymax),
                              n = 100, radius = 5/(2*Pi)){
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
    df[i,6:(6+length(tiles[[i]][[3]]))] <- tiles[[i]][[3]]
    df[i,17:(17+length(tiles[[i]][[4]]))] <- tiles[[i]][[4]]
  }
}
