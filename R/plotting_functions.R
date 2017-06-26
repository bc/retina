
##' @title Plot falciform process
##' @description smooths and plots the falciform process
##' @author Brian Cohn
##' @param falciform_x numeric vector of x coordinates
##' @param falciform_y numeric vector of y coordinates
plot_falciform_process <- function(falciform_x, falciform_y){
	fc_smoothed <- spline.poly(cbind(falciform_x,falciform_y), 50, k=10)
	polygon(fc_smoothed[,1], fc_smoothed[,2], col=rgb(0, 0, 0,0.5), lty="solid", border="gray42")
}


##' @title plot points of XY
##' @description visualize xy points
##' @author Brian Cohn
##' @param falciform_x numeric vector of x coordinates
##' @param falciform_y numeric vector of y coordinates
  plot_original_xy_locations <- function(x,y) {
  	points(x+0.005,y,pch=4, cex=0.5, col="gainsboro")
	points(x,y,pch=4, cex=0.5, col="black")
  }
##' @title Interpolate Input Data with Thin Plate Spline
##' @description Interpolation
##' @author Brian Cohn
##' @param minitics values referring to the circle
##' @param x input variable 1
##' @param y input variable 2
##' @param z response variable
##' @param lambda TPS parameter
##' @param polynomial_m TPS parameter
##' @param extrapolate true/false whether we should interpolate past the points, all the way to the eye equator


##' @title Plot the degree labels for latitudes
##' @description Plot degree numbers
##' @author Brian Cohn
##' @param outer.radius the extent of the radius of the plot
##' @param circle.rads the number of radian circles that are drawn
plot_degree_label_for_latitudes <- function(outer.radius, circle.rads) {
		axis(2, pos = -1.25 * outer.radius, at = sort(union(circle.rads,-circle.rads)), labels = NA)
		text( -1.26 * outer.radius, 
			sort(union(circle.rads, -circle.rads)),
			sort(union(circle.rads, -circle.rads)), 
			xpd = TRUE, pos = 2,family = "Palatino")
}



##' @title add a basic legend
##' @description plot a basic legend
##' @author Brian Cohn
##' @import fields
##' @param col Color vector
##' @param zlim limits of the densities
add_legend <- function(col, zlim) {
	par(mai = c(1,1,1.5,1.5))
	fields::image.plot(legend.only=TRUE, col=col, zlim=zlim, family = "Palatino")
}

##' @title draw line segments
##' @description put the radial lines on the plot
##' @author Brian Cohn
##' @param endpoints a 4 element numeric vector descrbing the xy and x'y' for the line segment.
##' @param color_hex hex string
draw_line_segments <- function(endpoints, color_hex="#66666650"){
  segments(endpoints[1], endpoints[2], endpoints[3], endpoints[4], col = color_hex )
}


##' @title Write labels at endpoint locations
##' @description put labels around the circle at each of the lines
##' @author Brian Cohn
##' @param r_label string of right label
##' @param l_label string of left label
##' @param degree the degree that is being placed in
##' @param endpoints vector of 4 numerics, x,y and x',y' defining the line segment
write_labels_at_endpoint_locations <- function(r_label, l_label, degree, endpoints){
		  lab1 <- bquote(.(r_label) * degree)
		  lab2 <- bquote(.(l_label) * degree)
		  text(endpoints[1], endpoints[2], lab1, xpd = TRUE, family = "Palatino")
		  text(endpoints[3], endpoints[4], lab2, xpd = TRUE, family = "Palatino")
	}


##' @title Remove points that are outside of the plotting circle
##' @description We do not need points plotted in the corners of the plotted circle.
##' @author Brian Cohn
##' @param minitics Spherical limit info
##' @param spatial_res spatial width in pixels of the plotted image
##' @param Mat Matrix of predicted points on the set grid
##' @param outer.radius max value of the radius
##' @param falciform_y numeric vector of y coordinates
nullify_vals_outside_the_circle <- function(minitics, spatial_res, heatmap_matrix, outer.radius){
  matrix_position_is_within_the_circle <- function() {!sqrt(markNA ^ 2 + t(markNA) ^ 2) < outer.radius}
  markNA <- matrix(minitics, ncol = spatial_res, nrow = spatial_res) 
  matrix_to_mask <- heatmap_matrix #matrix_to_mask is a mutable variable
  matrix_to_mask[matrix_position_is_within_the_circle()] <- NA #MUTABLE
  return(matrix_to_mask)
}


##' @title Compute Longitude Label Location
##' @description Find the location to put the label around the circle at each of the lines
##' @author Brian Cohn
##' @param axis.rads axis.rads
##' @param outer.radius numeric value for radius limit
##' @return label_locations the computed location for a label
compute_longitude_label_location <- function(axis.rads, outer.radius){
	return(c(RMat(axis.rads) %*% matrix(c(1.1, 0, -1.1, 0) * outer.radius, ncol = 2)))
}

##' @title Plot longitudinal spoke lines
##' @description put lines across the plotting circle
##' @author Brian Cohn
##' @param axis.radian radian
##' @param outer.radius numeric value of the outer radius limit
plot_longitudinal_lines <- function(axis_radian, outer.radius){
  endpoints <- zapsmall(c(RMat(axis_radian) %*% matrix(c(1, 0, -1, 0) * outer.radius,ncol = 2)))
  draw_line_segments(endpoints)
}

##' @title Plot longitudinal labels
##' @description put a degree label at the ends of the endpoints for each longitude
##' @author Brian Cohn
##' @param axis.rad axis radian
##' @param outer.radius numeric value of the outer radius limit
##' @param r_label label number
##' @param l_label label number
##' @param degree numeric, the degree of interest
plot_longitudinal_labels <- function(axis.rad, outer.radius, r_label, l_label, degree) {
  write_labels_at_endpoint_locations(r_label, l_label, degree, compute_longitude_label_location(axis.rad, outer.radius))
}


##' @title Internal equation for axis markup
##' @description define locations for radians with trigonometry
##' @author Brian Cohn
##' @param radians vector of radians
##' @return RMat trigonometric positions
RMat <- function(radians){
  return(matrix(c(cos(radians), sin(radians), -sin(radians), cos(radians)), ncol = 2))
}  


##' @title Draw Latitude Markings
##' @description plots radial lines, degree label for latitudes, and plots radial spokes with labels
##' @author Brian Cohn
##' @param radians vector of radians
##' @param circle.rads see fitplotazimuthal
##' @param outer.radius see fitplotazimuthal
draw_latitude_markings <- function(circle.rads, outer.radius) {
	if (missing(circle.rads)){circle.rads <- pretty(c(0,outer.radius-.4))}
	plot_circle_radial_lines(circle.rads)
	plot_degree_label_for_latitudes(outer.radius, circle.rads)
	plot_radial_spokes_and_labels(outer.radius)
}

##' @title Plot radial axes
##' @description Put radial axes onto visualization
##' @author Brian Cohn
##' @param outer.radius numeric value of the outer radius limit
plot_radial_spokes_and_labels <- function(outer.radius) {
	axis.rads <- c(0, pi / 6, pi / 3, pi / 2, 2 * pi / 3, 5 * pi / 6)
	r.labs <- c(90, 60, 30, 0, 330, 300)
	l.labs <- c(270, 240, 210, 180, 150, 120)
	for (i in 1:length(axis.rads)){ 
	  plot_longitudinal_lines(axis.rads[i], outer.radius)
	  plot_longitudinal_labels(axis.rads[i], outer.radius, r.labs[i], l.labs[i], degree)
	}
}


##' @title Polar Interpolation
##' @description This function will make a plot. Code from http://stackoverflow.com/questions/10856882/r-interpolated-polar-contour-plot was highly modified to meet retinal plotting funtionality.
##' @param x,y,z cartesian input in azimuthal format
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
##' @param outer.radius size of plot
##' @param circle.rads radius lines
##' @param spatial_res Used to define a spatial_res by spatial_res plotting resolution.
##' @param lambda lambda value for thin plate spline interpolation
##' @param xyrelief scaling factor for interpolation matrix.
##' @param tmp_input tmp_input
##' @param plot_suppress by default FALSE
##' @param compute_error whether to use fields::predictSE
##' @param falciform_coords vertices in xy format of the falciform process
##' @param falc2 a second falficorm coordinate file
##' @param polynomial_m A polynomial function of degree (m-1) will be included in the model as the drift (or spatial trend) component. Default is the value such that 2m-d is greater than zero where d is the dimension of x.
##' @param ... passed arguments
##' @import fields rgl RColorBrewer
##' @export
fit_plot_azimuthal<- function(
  ### Plotting data already post-azimuthal transformation
  x, y, z, 
  ### Plot component flags
  contours=TRUE,   # Add contours to the plotted surface
  legend=TRUE,        # Plot a surface data legend?
  axes=TRUE,      # Plot axes?
  should_plot_points=TRUE,        # Plot individual data points
  extrapolate=TRUE, # Should we extrapolate outside data points?
  ### Data splitting params for color scale and contours
  col_breaks_source = 2, # Where to calculate the color breaks from (1=data,2=surface)
  # If you know the levels, input directly (i.e. c(0,1))
  col_levels = 50,    # Number of color levels to use - must match length(col) if 
  #col specified separately
  col = rev(colorRampPalette(brewer.pal(11,"PuOr"))(col_levels)),  # Colors to plot
  contour_breaks_source = 1, # 1=z data, 2=calculated surface data
  # If you know the levels, input directly (i.e. c(0,1))<-default
  contour_levels = col_levels+1,
  ### Plotting params
  outer.radius = pi/2.0,  
  circle.rads = pretty(c(0,outer.radius)), #Radius lines
  spatial_res=1000, #Resolution of fitted surface
  single_point_overlay=0, #Overlay "key" data point with square 
  #(0 = No, Other = number of pt)
  ### Fitting parameters
  lambda=0.001, xyrelief=1,tmp_input = NULL, plot_suppress=FALSE,
  compute_error=FALSE,
  falciform_coords=NULL,
  falc2=NA, 
  polynomial_m=NULL,
  ...){ 

  minitics <- seq(-outer.radius, outer.radius, length.out = spatial_res)
  # interpolate the data

  vals <- interpolate_input_data(minitics, x, y, z, lambda, polynomial_m, extrapolate)
  t <- vals$t ; tmp <- vals$tmp; Mat <- vals$Mat;
  if (compute_error) {
  	error <- compute_thin_plate_spline_error(x,y,tmp)
  } else {
  	error <- NULL
  }

  if (plot_suppress == TRUE){
		return(list(t,tmp, error))
  }
  heatmap_matrix <- nullify_vals_outside_the_circle(minitics, spatial_res, Mat, outer.radius)
  
  zlim <- define_color_breaks_based_on_source(col_breaks_source,z, heatmap_matrix)

  init_square_mat_plot(heatmap_matrix, zlim, minitics, col)
  
  if (contours){ add_contours(minitics, heatmap_matrix,
  	contour_breaks=define_contour_breaks(contour_breaks_source, z, contour_levels, heatmap_matrix), xy)}
  
  plot_falciform_process(falciform_coords$x, falciform_coords$y)
  if (!is.na(falc2)) plot_falciform_process(falc2$x, falc2$y) #Plug in the secondary plot if it is available
  if (should_plot_points) plot_original_xy_locations(x,y)
  if (axes) draw_latitude_markings(circle.rads, outer.radius)
  if (legend) add_legend(col, zlim)

return(list(t,tmp, error))
}
