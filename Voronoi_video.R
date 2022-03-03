# Load libraries
library(ggvoronoi)
library(gganimate)
library(dplyr)

# Load data source and give columns names
ptshist_df<-trhistpts(histpts)


# Filtering for just first frame
ff_total <- filter(ptshist_df, Frame == 1)

# Defining pitch size for voronoi plot
rectangle <- data.frame(x=c(xmin,xmin,xmax+2*wid,xmax+2*wid),y=c(ymin,ymax,ymax,ymin))

# Make first frame and test what image looks like
ff <- ggplot(ff_total,aes(x,y)) +
  geom_voronoi(aes(fill=as.factor(y)),size=.125, outline = rectangle,show.legend = FALSE) +
  geom_vline(xintercept = xmax,color = 'white',linetype='solid',size=0.5) +
  geom_vline(xintercept = xmax+wid,color = 'white',linetype='solid',size=0.5) +
  stat_voronoi(geom="path",outline = rectangle) +
  geom_point(size=3) +
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
    ,legend.position="bottom"
  ) +
  labs(title = "      Evolution of the Voronoi Tesselation" # Title text
       ,subtitle = "         Metropolis' Algorithm"
       ,caption = "Author: Eloy Serrano       ")

# Create test image
ggsave(filename = paste0("frame_",1,".png") # filename
       ,plot = ff # variable for file
       ,width = 10, height = 4, dpi = 300, units = "in") # dimensions and image quality
#This will generate an image "test.png" which you can check whether this is the output you want for you gif, if not amend the code above.


# Loop through all frames
for(i in (min(ptshist_df$Frame):max(ptshist_df$Frame))){
  
  frame <- filter(ptshist_df, Frame == 1)
  
    plot <- ggplot(frame,aes(x,y)) +
    geom_voronoi(aes(fill=as.factor(y)),size=.125, outline = rectangle,show.legend = FALSE) +
    geom_vline(xintercept = xmax,color = 'white',linetype='solid',size=0.5) +
    geom_vline(xintercept = xmax+wid,color = 'white',linetype='solid',size=0.5) +
    stat_voronoi(geom="path",outline = rectangle) +
    geom_point(size=3) +
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
      ,legend.position="bottom"
    ) +
    labs(title = "      Evolution of the Voronoi Tesselation" # Title text
         ,subtitle = "         Metropolis' Algorithm"
         ,caption = "Author: Eloy Serrano       ")
  
  # Create test image
  ggsave(filename = paste0("frame_",i,".png") # filename
         ,plot = plot # variable for file
         ,width = 10, height = 4, dpi = 300, units = "in") # dimensions and image quality
  #This will generate an image "test.png" which you can check whether this is the output you want for you gif, if not amend the code above.
  
  # Writing print to console so you can check where the code is
  print(paste0("Frame: ",i," done!"))
}
