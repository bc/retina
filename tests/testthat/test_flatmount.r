context('flatmount for comparison')
test_that('flatmounts can be produced for comparison and inclusion in supplementary material', {
    
    gen_flatplot <- function(datapoints, xyz, outline_path, spatial_res = 1000){
	    colnames(datapoints) <- c("x", "y")
	    colnames(xyz) <- c("x_coord", "ycoord", "cells")
	    df <- cbind(datapoints[1:nrow(xyz),], xyz)
	    xy_image_coords <- RImageJROI::read.ijroi(outline_path)$coords
	    plot(df[,1],df[,2], xlab="x", ylab="y")
	    polygon(xy_image_coords[,1], xy_image_coords[,2])
	    fit <- Tps(cbind(df[,1], df[,2]), df$cells, lambda = 1e-5, m = 6)
	    bounds <- apply(xy_image_coords,2,range)
	    my_min <- 0
	    my_max <- max(bounds[2,])
	    minitics <- seq(my_min, my_max, length.out = spatial_res)

	    gridlist <- list(x = minitics, y = minitics) 
	    tmp <- predictSurface(fit, gridlist, extrap = FALSE)
	    library(raster)
	    rotate <- function(x) t(apply(x, 2, rev))
		r <- raster(rotate(tmp$z %>% reflect_across_horizontal_line))
		return(r)
    }

	# Use the full paths on your local computer
    datapoints <- read.csv("/Users/briancohn/Documents/GitHub/bc/retina/tests/testthat/test_retinas/raxp91/diagram_retina/datapoints.csv", header=TRUE)
    xyz <- read.csv("/Users/briancohn/Documents/GitHub/bc/retina/tests/testthat/test_retinas/raxp91/diagram_retina/xyz.csv", header=TRUE)
    outline_path <- "/Users/briancohn/Documents/GitHub/bc/retina/tests/testthat/test_retinas/raxp91/diagram_retina/outline.roi"
    
    surface <- gen_flatplot(datapoints, xyz, outline_path)
    plot(surface)
})