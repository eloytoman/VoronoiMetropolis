# Library rgl
library(rgl)

#Choose the size of the image on the output (800,650 to have 800 x 600)
r3dDefaults$windowRect <- c(0,50, 800, 650) 
open3d()
#If you want to put line on the background
#bg3d(sphere = TRUE, color = c("grey", "white"), lit = TRUE, back = "lines" ,lwd=2)

# This is to output a rgl plot in a rmarkdown document
# rgl::setupKnitr()


# plot
bg3d( col=rgb(0.2,0.8,0.5,0.8) )
theta <- seq(0, 2*pi, len = 50)
lon <- seq(0,20, len = 50)
zero <- seq(0,0, len = 50)
knot <- cylinder3d(
  center = cbind(zero,lon,zero),
  radius = 5/(2*pi),
  closed = FALSE)
shade3d(addNormals(subdivision3d(knot, depth = 2)), col = rgb(0.4,0.2,0.8,0.3))

# To display in an R Markdown document:
# rglwidget()

# save it as png
# snapshot3d( "~/Desktop/#20_portfolio_knot_3D.png", fmt="png")

# To save interactive plot to a file:
#htmlwidgets::saveWidget(rglwidget(width = 500, height = 500), 
#                        file = "HtmlWidget/3dknot.html",
#                        libdir = "libs",
#                        selfcontained = FALSE
#)


library(tessellation)
d <- delaunay(centricCuboctahedron())
v <- voronoi(d)
cell13 <- v[[13]]
isBoundedCell(cell13) # TRUE
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
invisible(lapply(cell13[["cell"]], function(edge){
edge$plot(edgeAsTube = TRUE, tubeRadius = 0.025, tubeColor = "yellow")
}))
cellvertices <- cellVertices(cell13)
spheres3d(cellvertices, radius = 0.1, color = "green")


library(tessellation)
d <- delaunay(centricCuboctahedron())
v <- voronoi(d)
cell13 <- v[[13]]
isBoundedCell(cell13) # TRUE
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
plotBoundedCell3D(
  cell13, edgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "yellow",
  facetsColor = "navy", alpha = 0.7
)



library(tessellation)
points <- rbind(
  c(0.5,0.5,0.5),
  c(0,0,0),
  c(0,0,1),
  c(0,1,0),
  c(0,1,1),
  c(1,0,0),
  c(1,0,1),
  c(1,1,0),
  c(1,1,1)
)
#points<-cbind(Radius*cos(x1),Radius*sin(x1),y1)
del <- delaunay(points)
del$vertices[[1]]
del$tiles[[1]]
del$tilefacets[[1]]
# an elevated Delaunay tessellation ####
f <- function(x, y){
  dnorm(x) * dnorm(y)
}
x <- y <- seq(-5, 5, length.out = 50)
grd <- expand.grid(x = x, y = y) # grid on the xy-plane
points <- as.matrix(transform( # data (x_i, y_i, z_i)
  grd, z = f(x, y)
))
del <- delaunay(points, elevation = TRUE)
del[["volume"]] # close to 1, as expected
# plotting
library(rgl)
mesh <- del[["mesh"]]
open3d(windowRect = c(100, 100, 612, 356), zoom = 0.6)
aspect3d(1, 1, 20)
shade3d(mesh, color = "limegreen")
wire3d(mesh)


tetrahedron <-
  rbind(
    c(2*sqrt(2)/3, 0, -1/3),
    c(-sqrt(2)/3, sqrt(2/3), -1/3),
    c(-sqrt(2)/3, -sqrt(2/3), -1/3),
    c(0, 0, 1)
  )
angles <- seq(0, 2*pi, length.out = 91)[-1]
R <- 2.5
circle1 <- t(vapply(angles, function(a) R*c(cos(a), sin(a), 0), numeric(3L)))
circle2 <- t(vapply(angles, function(a) R*c(cos(a), 0, sin(a)), numeric(3L)))
circle3 <- t(vapply(angles, function(a) R*c(0, cos(a), sin(a)), numeric(3L)))
circles <- rbind(circle1, circle2, circle3)
pts <- rbind(tetrahedron, circles)
points<-cbind(Radius*cos(x1),Radius*sin(x1),y1)
d <- delaunay(points, degenerate = TRUE)
v <- voronoi(d)
open3d(windowRect = c(50, 50, 562, 562))
material3d(lwd = 2)
plotVoronoiDiagram(v, luminosity = "bright")
