
##' Splining a polygon
##' @description
# .  The rows of 'xy' give coordinates of the boundary
# vertices, in order.  'vertices' is the number of spline vertices to create.
# (Not all are used: some are clipped from the ends.)  'k' is the number of
# points to wrap around the ends to obtain a smooth periodic spline.  reference:
# http://gis.stackexchange.com/questions/24827/how-to-smooth-the-polygons-in-a-contour-map
# Returns an array of points.
##' @param xy xy dataframe coordinates, in two columns.
##' @param k number of points to wrap around ends to obtain smooth spline
##' @param ... further arguments to spline
##' @return xy xy dataframe coordinates, in two columns, with new points from spline.
##' @author Brian Cohn \email{brian.cohn@@usc.edu}, Lars Schmitz
##' @export
##' @importFrom stats spline
spline_poly <- function(xy, vertices, k = 3, ...) {
    # Assert: xy is an n by 2 matrix with n >= k.

    # Wrap k vertices around each end.
    n <- dim(xy)[1]
    if (k >= 1) {
        data <- rbind(xy[(n - k + 1):n, ], xy, xy[1:k, ])
    } else {
        data <- xy
    }

    # Spline the x and y coordinates.
    data.spline <- spline(1:(n + 2 * k), data[, 1], n = vertices, ...)
    x <- data.spline$x
    x1 <- data.spline$y
    x2 <- spline(1:(n + 2 * k), data[, 2], n = vertices, ...)$y

    # Retain only the middle part.
    cbind(x1, x2)[k < x & x <= n + k, ]
}
