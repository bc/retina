require(rgl)
require(fields)
require(RColorBrewer)
require(sphereplot)
require(mapproj)

##' @title Plot from matrix
##' @description Sometimes users will want to permute their density matrices prior to plotting. This can be useful when multiplying the densities by some scaling factor, accommodating for bias, or incorporating additional noise for analysis with robustness.
##' @param contours whether to plot contours.
##' @param legend Color legend with tick marks
##' @param axes Radial axes
##' @param points whether to plot individual datapoints as X's
##' @param extrapolate By default TRUE, will make plot a circle.
##' @param col_breaks_source 2 element vector with max and min
##' @param col_levels number of color levels
##' @param col colors to plot
##' @param contour_breaks_source 1 if data, 2 if calculated surface data
##' @param contour_levels number of contour levels
##' @param outer_radius size of plot
##' @param circle.rads radius lines
##' @param spatial_res Used to define a spatial_res by spatial_res plotting resolution.
##' @param single_point_overlay Overlay "key" data point with square (0 = No, Other = number of pt)
##' @param interp.type depreciated
##' @param lambda lambda value for thin plate spline interpolation
##' @param xyrelief scaling factor for interpolation matrix.
##' @param z1 Vector of density values (used to get the max and min)
##' @param z2 Vector of density values (used to get the max and min)
##' @param MATRIX Matrix object containing density values at each row and col
##' @param ... further arguments passed to or from other methods.
##' @details Takes in a square matrix of azimuthal equidistant plot projection points. Generates plot projection.
##' @export
plot_from_MAT <- function(
  z1, z2, MATRIX,
  contours=TRUE,
  legend=TRUE,
  axes=TRUE,
  points=TRUE,
  extrapolate=FALSE,
  col_breaks_source = 2,
  col_levels = 10,
  col = rev(colorRampPalette(brewer.pal(11,"PuOr"))(col_levels)),
  contour_breaks_source = 1,
  contour_levels = col_levels+1,
  outer_radius = pi/2,
  circle.rads = pretty(c(0,outer_radius)),
  spatial_res=128,
  single_point_overlay=0,
  interp.type = 1,
  lambda=0.001, xyrelief=1,...) {
  z <- c(min(MATRIX), max(MATRIX))
  spatial_res <- length(MATRIX[,1])
# interpolate the data


  if (interp.type == 1){

	minitics <- seq(-outer_radius, outer_radius, length.out = spatial_res)
	grid.list = list(x=minitics,y=minitics) #choose locations to predict at
	Mat = MATRIX
  }
  else {stop("interp.type value not valid")}





#   ######START
#   ###;; ________________________________________________________________________
  # mark cells outside circle as NA
  markNA <- matrix(minitics, ncol = spatial_res, nrow = spatial_res)
  Mat[!sqrt(markNA ^ 2 + t(markNA) ^ 2) < outer_radius] <- NA

  ### Set contour_breaks based on requested source
  if ((length(contour_breaks_source == 1)) & (contour_breaks_source[1] == 1)){
	contour_breaks = seq(min(z,na.rm=TRUE),max(z,na.rm=TRUE),
						 by=(max(z,na.rm=TRUE)-min(z,na.rm=TRUE))/(contour_levels-1))
  }
  else if ((length(contour_breaks_source == 1)) & (contour_breaks_source[1] == 2)){
	contour_breaks = seq(min(Mat,na.rm=TRUE),max(Mat,na.rm=TRUE),
						 by=(max(Mat,na.rm=TRUE)-min(Mat,na.rm=TRUE))/(contour_levels-1))
  }
  else if ((length(contour_breaks_source) == 2) & (is.numeric(contour_breaks_source))){
	print(paste0("Manual contour range set from ", contour_breaks_source[1], " to ", contour_breaks_source[2]))
	contour_breaks = pretty(contour_breaks_source,n=contour_levels)
	contour_breaks = seq(contour_breaks_source[1],contour_breaks_source[2],
						 by=(contour_breaks_source[2]-contour_breaks_source[1])/(contour_levels-1))
  }
  else {stop("Invalid selection for \"contour_breaks_source\"")}

  ### Set color breaks based on requested source
  if ((length(col_breaks_source) == 1) & (col_breaks_source[1] == 1))
  {zlim=c(min(z,na.rm=TRUE),max(z,na.rm=TRUE))}
  else if ((length(col_breaks_source) == 1) & (col_breaks_source[1] == 2))
  {zlim=c(min(Mat,na.rm=TRUE),max(Mat,na.rm=TRUE))}
  else if ((length(col_breaks_source) == 2) & (is.numeric(col_breaks_source)))
  {zlim=col_breaks_source}
  else {stop("Invalid selection for \"col_breaks_source\"")}


  # begin plot
  Mat_plot = Mat
  Mat_plot[which(Mat_plot<zlim[1])]=zlim[1]
  Mat_plot[which(Mat_plot>zlim[2])]=zlim[2]
  image(x = minitics, y = minitics, Mat_plot , useRaster = TRUE, asp = 1, axes = FALSE, xlab = "", ylab = "", zlim = zlim, col = col)

  # add contours if desired
  if (contours){
	CL <- contourLines(x = minitics, y = minitics, Mat, levels = contour_breaks)
	A <- lapply(CL, function(xy){
	  lines(xy$x, xy$y, col = gray(.2), lwd = .5)
	})
  }

  # add radial axes if desired
  if (axes){
	# internals for axis markup
	RMat <- function(radians){
	  matrix(c(cos(radians), sin(radians), -sin(radians), cos(radians)), ncol = 2)
	}

	circle <- function(x, y, rad = 1, nvert = 500){
	  rads <- seq(0,2*pi,length.out = nvert)
	  xcoords <- cos(rads) * rad + x
	  ycoords <- sin(rads) * rad + y
	  cbind(xcoords, ycoords)
	}

	# draw circles
	if (missing(circle.rads)){
	  circle.rads <- pretty(c(0,outer_radius-.4))
	}

	for (i in circle.rads){
	  lines(circle(0, 0, i), col = "#66666650")
	}

	# put on radial spoke axes:
	axis.rads <- c(0, pi / 6, pi / 3, pi / 2, 2 * pi / 3, 5 * pi / 6)
	r.labs <- c(90, 60, 30, 0, 330, 300)
	l.labs <- c(270, 240, 210, 180, 150, 120)

	for (i in 1:length(axis.rads)){
	  endpoints <- zapsmall(c(RMat(axis.rads[i]) %*% matrix(c(1, 0, -1, 0) * outer_radius,ncol = 2)))
	  segments(endpoints[1], endpoints[2], endpoints[3], endpoints[4], col = "#66666650")
	  endpoints <- c(RMat(axis.rads[i]) %*% matrix(c(1.1, 0, -1.1, 0) * outer_radius, ncol = 2))
	  lab1 <- bquote(.(r.labs[i]) * degree)
	  lab2 <- bquote(.(l.labs[i]) * degree)
	  text(endpoints[1], endpoints[2], lab1, xpd = TRUE)
	  text(endpoints[3], endpoints[4], lab2, xpd = TRUE)
	}
	axis(2, pos = -1.25 * outer_radius, at = sort(union(circle.rads,-circle.rads)), labels = NA)
	text( -1.26 * outer_radius, sort(union(circle.rads, -circle.rads)),sort(union(circle.rads, -circle.rads)), xpd = TRUE, pos = 2)
  }
  if (legend){

	image.plot(legend.only=TRUE, col=col, zlim=zlim, smallplot=c(0.81,.84,0.2,.8))
  }
}
