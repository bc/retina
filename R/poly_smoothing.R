
##' Splining a polygon
##' @description
##' The rows of 'xy_dataframe' give coordinates of the boundary
##' vertices, in order.  'vertices' is the number of spline vertices to create.
##' (Not all are used: some are clipped from the ends.)  'k' is the number of
##' points to wrap around the ends to obtain a smooth periodic spline.  reference:
##' http://gis.stackexchange.com/questions/24827/how-to-smooth-the-polygons-in-a-contour-map
##' Returns an array of points.
##' @param xy_dataframe xy dataframe coordinates, in two columns.
##' @param k number of points to wrap around ends to obtain smooth spline
##' @param ... further arguments to spline
##' @return xy_dataframe_splined xy dataframe coordinates, in two columns, with new points from spline.
##' @author Brian Cohn \email{brian.cohn@@usc.edu}, Lars Schmitz
##' @export
##' @importFrom stats spline
spline_poly <- function(xy_dataframe, vertices, k = 3, ...) {
    # Wrap k vertices around each end.
    n <- dim(xy_dataframe)[1]
    if (k >= 1) {
        data <- rbind(xy_dataframe[(n - k + 1):n, ], xy_dataframe, xy_dataframe[1:k, ])
    } else {
        data <- xy_dataframe
    }
    # Spline the x and y coordinates.
    data.spline <- spline(1:(n + 2 * k), data[, 1], n = vertices, ...)
    x <- data.spline$x
    x1 <- data.spline$y
    x2 <- spline(1:(n + 2 * k), data[, 2], n = vertices, ...)$y

    # Retain only the middle part.
    return(cbind(x1, x2)[k < x & x <= n + k, ])
}
