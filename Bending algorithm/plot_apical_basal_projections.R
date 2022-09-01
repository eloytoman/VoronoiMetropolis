library(deldir)
library(ggplot2)
library(ggvoronoi)
library(dplyr)

ggplot_vororonoi_analysis <-function(pointsapical, pointsbasal, ratio = 2.5){
  
  rectangle <- data.frame(x = c(xmin, xmin, xmax+2*wid, xmax+2*wid),
                          y = c(ymin, ymax, ymax, ymin))
  rectanglebasal <- data.frame(x = c(xmin, xmin, ratio*(xmax + 2*wid), ratio*(xmax+2*wid)), 
                               y = c(ymin, ymax, ymax, ymin))
  pl <- ggplot(pointsapical, aes(x, y)) +
    geom_voronoi(aes(fill = as.factor(pt)), size=.125, outline = rectangle, show.legend = FALSE) +
    geom_vline(xintercept = xmax,color = 'white', linetype='solid', size=1) +
    geom_vline(xintercept = xmax+wid, color = 'white', linetype='solid', size=1) +
    stat_voronoi(geom="path", outline = rectangle) +
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
    labs(title = "Apical surface tessellation" # Title text
         ,caption = "Author: Eloy Serrano        ")
  show(pl)
  
  pl2 <- ggplot(pointsbasal,aes(x,y)) +
    geom_voronoi(aes(fill = as.factor(pt)), size=.125, outline = rectanglebasal, show.legend = FALSE) +
    geom_vline(xintercept = xmax*ratio,color = 'white',linetype='solid',size=1) +
    geom_vline(xintercept = (xmax+wid)*ratio,color = 'white',linetype='solid',size=1) +
    stat_voronoi(geom="path",outline = rectanglebasal) +
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
    labs(title = "Basal surface tessellation" # Title text
         ,caption = "Author: Eloy Serrano        ")
  show(pl2)
  
  pointsapicalproy <- pointsapical
  pointsapicalproy$x <- pointsapicalproy$x*ratio
  
  pl3 <- ggplot(pointsapicalproy,aes(x,y)) +
    geom_voronoi(aes(fill = as.factor(pt)), size=.125, outline = rectanglebasal, show.legend = FALSE) +
    geom_vline(xintercept = xmax*ratio,color = 'white',linetype='solid',size=1) +
    geom_vline(xintercept = (xmax+wid)*ratio,color = 'white',linetype='solid',size=1) +
    stat_voronoi(geom="path",outline = rectanglebasal) +
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
    labs(title = "Proyection of the apical surface tessellation" # Title text
         ,caption = "Author: Eloy Serrano        ")
  show(pl3)
}

xmin <- 0
xmax <- 5
ymin <- 0
ymax <- 20
ratio <- 2.5
wid <-5
histpts0<-histpts
pointsapa0 <- dplyr::filter(histpts0[[1]], Frame == 100 );
pointsbasa0 <- dplyr::filter(histpts0[[4]], Frame == 100 );
pointsapa0 <- pointsapa0[,-3]
pointsbasa0 <- pointsbasa0[,-3]
pointsapa0[,3] <- 1:100
pointsbasa0[,3]<- 1:100
names(pointsapa0) <- c("x","y","pt")
names(pointsbasa0) <- c("x","y","pt")

ggplot_vororonoi_analysis(pointsapical = points[[1]], pointsbasal = points[[12]])

# pointsapa5 <- dplyr::filter(histpts5[[1]], Frame == 100 );
# pointsbasa5 <- dplyr::filter(histpts5[[4]], Frame == 100 );
# pointsapa5 <- pointsapa5[,-3]
# pointsbasa5 <- pointsbasa5[,-3]
# pointsapa5[,3] <- 1:100
# pointsbasa5[,3]<- 1:100
# names(pointsapa5) <- c("x","y","pt")
# names(pointsbasa5) <- c("x","y","pt")






