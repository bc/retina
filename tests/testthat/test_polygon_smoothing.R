context("polygon (falciform) smoothing")

test_that("polygon is smoothed appropriately", {

#
# Example polygon, randomly generated.
#
set.seed(17)
n.vertices <- 10
theta <- (runif(n.vertices) + 1:n.vertices - 1) * 2 * pi / n.vertices
r <- rgamma(n.vertices, shape=3)
xy <- cbind(cos(theta) * r, sin(theta) * r)
s <- spline.poly(xy, 100, k=3)
plot(s, type="l", lwd=2, col="Gray")
lines(xy, col="Red", lty=2, lwd=2)
points(xy, col="Red", pch=19)
points(s, cex=0.8)
points(s[c(1,dim(s)[1]),], col="Blue", pch=19)

})