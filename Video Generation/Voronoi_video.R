# Load libraries
library(ggplot2)
library(ggvoronoi)
library(gganimate)
library(dplyr)
library(deldir)



byareaenergy<-function(points){
  tesel<-deldir(points$x,points$y,rw=rec)
  tilest<-tile.list(tesel)[(n+1):(2*n)]
  areas<-sapply(tilest,function(x){x$area})
  points[,4]<-c(areas,areas,areas)
  perim_ad<-(tilePerim(tilest)$perimeters)/sqrt(A0)
  areas_ad<-(areas)/A0
  encells<-(areas_ad-1)^2+(gam_ad/2)*perim_ad+lambda_ad*perim_ad
  points[,5]<-c(encells,encells,encells)
  
  names(points)<-c("x","y","Frame","area","energy")
  return(points)
}

ptshist_df<-data.frame(x=c(),y=c(),Frame=c(), area=c(),energy=c())

for (i in 1:100) {
  a<-length(ptshist_df$x)
  ptshist_df[(a+1):(a+3*n),c(1,2,3,4,5)]<-byareaenergy(dplyr::filter(histpts, Frame == i))
}

# Load data source and give columns names
names(ptshist_df)<-c("x","y","Frame","Area of the cell","Relative energy of the cell")

minar<-min(df$`Area of the cell`)
maxar<-max(df$`Area of the cell`)

# Filtering for just first frame
ff_total1 <- byareaenergy(filter(histpts, Frame == 1))
ff_total100 <- byareaenergy(filter(histpts, Frame == 100))

# Defining pitch size for voronoi plot
rectangle <- data.frame(x=c(xmin,xmin,xmax+2*wid,xmax+2*wid),y=c(ymin,ymax,ymax,ymin))

# Make first frame and test what image looks like
ff <- ggplot(ff_total1,aes(x,y)) +
  geom_voronoi(aes(fill=area), size=.125, outline = rectangle) +
  geom_vline(xintercept = xmax ,color = 'white',linetype='solid',size=0.5) +
  geom_vline(xintercept = xmax+wid,color = 'white',linetype='solid',size=0.5) +
  stat_voronoi(geom="path", outline = rectangle) +
  geom_point(size=3) +
  scale_fill_gradient(low = "dark blue", high = "white",
                      limits=c(min(ff_total1$area),
                               max(ff_total1$area)))+
  theme(
    panel.grid.major = element_blank() # Remove gridlines (major)
    ,panel.grid.minor = element_blank() # Remove gridlines (minor)
    ,panel.background = element_blank() # Remove grey background
    ,plot.title = element_text(hjust = 0, size = 20, colour = "#323232") # Title size and colour
    ,plot.subtitle = element_text(hjust = 0, size = 14, colour = "#323232") # Subtitle size and colour
    ,plot.caption = element_text(vjust = 0.3, size = 11, colour = "#323232") # Caption size and colour
    ,axis.ticks.y = element_blank() # Remove tick marks (Y-Axis)
    ,axis.text.y =  element_blank() # Remove scale marks (Y-Axis)
    ,axis.title.y = element_blank() # Remove axis label (Y-Axis) 
    ,axis.ticks.x = element_blank() # Remove tick marks (X-Axis)
    ,axis.text.x  = element_blank() # Remove axis scale (X-Axis)
    ,axis.title.x = element_blank() # Remove axis label (X-Axis) 
    ,legend.position="right"
  ) +
  labs(title = "      Evolution of the Voronoi Tesselation" # Title text
       ,subtitle = "         Metropolis' Algorithm"
       ,caption = "Author: Eloy Serrano       ")

show(ff)

# Create test image
ggsave(filename = paste0("frame_",1,".png") # filename
       ,plot = ff # variable for file
       ,width = 10, height = 7, dpi = 300, units = "in") # dimensions and image quality
#This will generate an image "test.png" which you can check whether this is the output you want for you gif, if not amend the code above.
save(filename = paste0("frame_",1,".png") # filename
     ,plot = ff # variable for file
     ,width = 10, height = 7, dpi = 300, units = "in") # dimensions and image quality
#This will generate an image "test.png" which you can check whether this is the output you want for you gif, if not amend the code above.



minframe<-min(ptshist_df$Frame)
# Loop through all frames
for(i in (min(ptshist_df$Frame):max(ptshist_df$Frame))){
  
  frame <- filter(ptshist_df, Frame == i)
  
    plot <- ggplot(frame,aes(x,y)) +
      geom_voronoi(aes(fill=`Area of the cell`),size=.125, outline = rectangle) +
      geom_vline(xintercept = xmax,color = 'white',linetype='solid',size=0.5) +
      geom_vline(xintercept = xmax+wid,color = 'white',linetype='solid',size=0.5) +
      stat_voronoi(geom="path",outline = rectangle) +
      geom_point(size=3) +
      scale_fill_gradient(low = "dark orange", high = "white", limits=c(minar,maxar))+
      theme(
        panel.grid.major = element_blank() # Remove gridlines (major)
        ,panel.grid.minor = element_blank() # Remove gridlines (minor)
        ,panel.background = element_blank() # Remove grey background
        ,plot.title = element_text(hjust = 0, size = 20, colour = "#323232") # Title size and colour
        ,plot.subtitle = element_text(hjust = 0, size = 14, colour = "#323232") # Subtitle size and colour
        ,plot.caption = element_text(vjust = 0.3, size = 11, colour = "#323232") # Caption size and colour
        ,axis.ticks.y = element_blank() # Remove tick marks (Y-Axis)
        ,axis.text.y =  element_blank() # Remove scale marks (Y-Axis)
        ,axis.title.y = element_blank() # Remove axis label (Y-Axis) 
        ,axis.ticks.x = element_blank() # Remove tick marks (X-Axis)
        ,axis.text.x  = element_blank() # Remove axis scale (X-Axis)
        ,axis.title.x = element_blank() # Remove axis label (X-Axis) 
        ,legend.position="right"
      ) +
      labs(title = "      Evolution of the Voronoi Tesselation" # Title text
           ,subtitle = "         Metropolis' Algorithm"
           ,caption = "Author: Eloy Serrano       ")
    
  # Create test image
  ggsave(filename = paste0("frame_",i,".png") # filename
         ,plot = plot # variable for file
         ,width = 10, height = 7, dpi = 300, units = "in") # dimensions and image quality
  #This will generate an image "test.png" which you can check whether this is the output you want for you gif, if not amend the code above.
  
  # Writing print to console so you can check where the code is
  print(paste0("Frame: ",i," done!"))
}
