
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
##' @param radius_vals vector of radius values (for latitudes)
##' @param outer.radius see fitplotazimuthal
draw_latitude_markings <- function(radius_vals, outer.radius) {
	plot_circle_radial_lines(radius_vals)
	plot_degree_label_for_latitudes(outer.radius, radius_vals)
	plot_radial_spokes_and_labels(outer.radius)
}


##' @title Draw circle radians about the center
##' @description Draw N radius lines (circles)
##' @author Brian Cohn
##' @param radius_vals vector of radius values where the circles will be drawn
plot_circle_radial_lines <- function(radius_vals, color_hex = "#66666650"){
	circle <- function(x, y, rad = 1, nvert = 500){
	  rads <- seq(0,2*pi,length.out = nvert)
	  xcoords <- cos(rads) * rad + x
	  ycoords <- sin(rads) * rad + y
	  cbind(xcoords, ycoords)
	}
	for (i in radius_vals){
	  lines(circle(0, 0, i), col = color_hex)
	}
}

##' @title Define Zlim by the requested color breaks source
##' @description Set color breaks (zlim) based on requested source
##' @author Brian Cohn
##' @param col_breaks_source A 2 element vector with max and min
##' @param z the response values from the input data (the retinal densities)
##' @param Mat The predicted retinal densities across the xy space
define_color_breaks_based_on_source <- function(col_breaks_source,z, Mat) {
  if ((length(col_breaks_source) == 1) & (col_breaks_source[1] == 1))
  	{zlim <- c(min(z,na.rm=TRUE),max(z,na.rm=TRUE))}
  else if ((length(col_breaks_source) == 1) & (col_breaks_source[1] == 2))
  	{zlim <- c(min(Mat,na.rm=TRUE),max(Mat,na.rm=TRUE))}
  else if ((length(col_breaks_source) == 2) & (is.numeric(col_breaks_source)))
  	{zlim <- col_breaks_source}
  else {stop("Invalid selection for \"col_breaks_source\"")}
  return(zlim)
}

 interpolate_input_data <- function(minitics, x, y, z, lambda, polynomial_m, extrapolate){
	grid.list = list(x=minitics,y=minitics) #choose locations to predict at
	t <- Tps(cbind(x,y),z,lambda=lambda, m = polynomial_m) #computationally intensive
	tmp <- predictSurface(t,grid.list,extrap=extrapolate)
	Mat = tmp$z
	return(list(t=t, tmp=tmp, Mat= Mat))
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

##' @title Create a pretty number list up to but not including the upper limit
##' @description Remove the last element if it exceeds the upper limit
##' @author Brian Cohn
##' @param upper_limit numeric upper bound
##' @param lower_limit numeric lower bound
##' @param pretty_vec numeric vector
pretty_list_not_including_max <- function(lower_limit, upper_limit){
	radian_list <- pretty(c(lower_limit,upper_limit))
	if (max(radian_list) > upper_limit){
		return(a[1:length(a)-1])
	} else {
		return(radian_list)
	}
}

##' @title Polygon Spline Fit
##' @details enhance the resolution of a polygon verticies dataframe by creating a spline along each vertex.
##' @param xy vertices in dataframe with x and y columns, in order (not all are used).
##' @param vertices Number of spline vertices to create.
##' @param k Wraps K vertices around each end. n >=k
##' @param ... further arguments passed to or from other methods.
##' @return Coords More finely placed vertices for the polygon.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @references http://gis.stackexchange.com/questions/24827/how-to-smooth-the-polygons-in-a-contour-map
spline.poly <- function(xy, vertices, k=3, ...) {
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y

  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

##' @title Retina Plot
##'
##' @description
##' \code{retinaplot} Generates an Azimuthal Equidistant plot projection of a retina object.
##' You can also set lambda(floating point number) and polynomial_m(integer), as well as extrapolate (TRUE, FALSE).
##' @param inner_eye_view boolean, default is TRUE. If set to false, the plotted view of the retina will have the viewpoint within the skull looking at the rear of the eye. inner_eye_view has the same view as the traditional wholemount.
##' @param ... further arguments passed to or from other methods.
##' @param rotation degrees to rotate CCW (int or floating point)
##' @param return_fit logical, whether or not to return the interpolation fit data.
##' @param spatial_res define the number of pixels (resolution) the plot will be
##' @param retina_object A list containing an element \code{azimuthal_data.datapoints} with
##' \code{x,y,z} datapoints. File must also include \code{azimuthal_data.falciform}.
##' @return Base-R polar plot
##' 
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' 
##' 
##' @family visualization
##' @export
retinaplot <- function(retina_object, spatial_res=1000, rotation=0, inner_eye_view=TRUE, return_fit=FALSE, ...){
  AZx   <- retina_object$azimuthal_data.datapoints[[1]]$x
  AZy   <- retina_object$azimuthal_data.datapoints[[1]]$y
  AZz   <- retina_object$azimuthal_data.datapoints[[1]]$z
  # if (rotation !=0){
  rotAZ <- cartesian_rotation(AZx, AZy, rotation)
  AZx   <- rotAZ$x
  AZy   <- rotAZ$y
  rotFALC <- cartesian_rotation(retina_object$azimuthal_data.falciform[[1]]$x,
                  retina_object$azimuthal_data.falciform[[1]]$y, rotation)
  retina_object$azimuthal_data.falciform[[1]]$x <- rotFALC$x
  retina_object$azimuthal_data.falciform[[1]]$y <- rotFALC$y
  message(paste("rotated by", rotation, "degrees"))
  # }
  if (inner_eye_view==TRUE){
    AZx = AZx*-1.0
    retina_object$azimuthal_data.falciform[[1]]$x <- retina_object$azimuthal_data.falciform[[1]]$x*-1.0
  }
  temp <- fit_plot_azimuthal(
      AZx,
      AZy,
      AZz,
      outer.radius=1.6,
      spatial_res=spatial_res,
      falciform_coords=retina_object$azimuthal_data.falciform[[1]],
      falc2=NA, ...
      )
  if (return_fit){return(temp)}
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
  circle.rads = pretty_list_not_including_max(0,outer.radius), #Radius lines
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
