context('flatmount for comparison')
test_that('flatmounts can be produced for comparison and inclusion in supplementary material', {
    datapoints <- read.csv("/Users/briancohn/Documents/GitHub/bc/retina/tests/testthat/test_retinas/3hgbqg/diagram_retina/datapoints.csv", header=TRUE)
    colnames(datapoints) <- c("x", "y")
    xyz <- read.csv("/Users/briancohn/Documents/GitHub/bc/retina/tests/testthat/test_retinas/3hgbqg/diagram_retina/xyz.csv", header=TRUE)
    colnames(xyz) <- c("x_coord", "ycoord", "cells")
    df <- cbind(datapoints[1:nrow(xyz),], xyz)
    out <- Tps(expand.grid(df$x, df$y), df$z)
    ggplot(df, aes(x,y,z=cells))  + geom_raster(aes(fill = cells)) + geom_contour(colour = "white")
}