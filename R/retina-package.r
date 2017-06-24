##' retina.
##' @description
##' Copyright (c) 2013-2014 Brian Cohn and Lars Schmitz
##' @name retina
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @docType package
NULL

require(retistruct)
require(rgl)
require(fields)
require(RColorBrewer)
require(sphereplot)
require(mapproj)

##' Reef Fish
##'
##' A dataset containing three retinal objects, computed using the retina package: Pmol 753,Pmol 752 and Ntae 381.
##'
##'
##' @docType data
##' @keywords datasets
##' @name reef_fish
##' @usage data(reef_fish)
##' @format A list (in retina object form)
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
NULL

#' Pseudodax 753 Retinal object
#' @docType data
#' @keywords datasets
#' @format retina_object
#' @name Pmol_753
NULL

#' Pseudodax 752 Retinal object
#' @docType data
#' @keywords datasets
#' @format retina_object
#' @name Pmol_752
NULL

#' Ntae_381 Retinal object
#' @docType data
#' @keywords datasets
#' @format retina_object
#' @name Ntae_381	
NULL

##' @title Half ellipse perimeter approximation (ellipsoid assumption)
##' @description
##' Ramanujan approximation divided by two.
##' \eqn{p \approx \pi * [3(a+b)-\sqrt{(3a + b)(a + 3b)}]} in half
##' @param a Major axis (longer)
##' @param b Minor axis (shorter)
##' @return semielliptical perimeter in same units as input
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @export
semi_ellipse_perimeter <- function( a, b)
{
	p <- pi * (3*(a+b)- sqrt((3*a + b)*(a + 3*b)))
	return(p/2)
}

##' Retinal Perimeter Estimation (spherical assumption)
##' @description
##' Calculates the retinal perimeter distance in same units as input, with the assumption that the retinal shape is hemispherical.
##' @details Assumes a spherical eye, where eye diameter approximates twice the sphere's radius.
##' @param ED Eye Diameter of the eye sample, measured at the widest point of the eye.
##' @return Retinal Arclength in millimeters
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @export
retinal_arclen <- function(ED){
	C_eye <- ED*pi #assuming spherical eye with ED = sphere diameter
	retinal_arclen <- C_eye/2
	return(retinal_arclen)
}

##' Rim Latitude Estimation (Ellipsoidal assumption)
##' @details Assumes an ellipsoidal eye
##' @param ED Eye Diameter of the eye sample, measured at the widest points of the eye.
##' @param AL Axial Length of the eye sample, measured at the longest measure of the eye from anterior to posterior.
##' @return Retinal latitude in degrees
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @export
retinal_phi0 <- function(ED, AL){
	return(90-180/pi*acos(2*(AL/ED)-1)) ##This line written by manuscript reviewer via Journal of Vision
}

##' @title import_xyz
##' @description
##' Grabs the XYZ from the path
##' @param path Path to a folder containing a three column comma separated value database. File must be named \code{xyz.csv}
##' @return dataframe with X Y and Z columns with NA's. Each row is a density measurement at a given XY location in ImageJ coordinates.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family munge
import_xyz <- function(path){
	#Input path to folder containing xyz.csv, Reads and condenses count data per sampling site. Make sure you don't have a slash at the end
	xyz <- na.omit(read.csv(paste0(path,"/xyz.csv")))
	return(xyz)
}

##' Polynomial vs Lambda Visualization
##'
##' @description
##' Visualize the effect of polynomial number and lambda on the thin plate spline fit.
##' @param retina_obj Retinal object
##' @param spatial_res number of pixels wide the output image will be. (1000 by default)
##' @param frequency_cap y axis limit for the frequency (used to set all plots to the same Y scale)
##' @param polynomial_m_vec by default c(3,4), defines the polynomial degrees to try
##' @param lambda_vec by default c(0,0.001), defines the lambda values to try
##' @param error_x_limits xlim 2 element vector passed for the x axis. See ?plot
##' @details pdf documents in the working directory for each combination of lambda and m (polynomial degree)
##' @return DAT data.frame containing each column with a combination of lambda and polynomial values. Each row value is an error value
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu} Lars Schmitz
##' @export

polynomial_vs_lambda<- function(retina_obj, spatial_res= 1000, polynomial_m_vec=c(3,4), lambda_vec=c(0,0.001), frequency_cap=100, error_x_limits= c(1000,7000)){
	point_num <- length(retina_obj$fit_data1$y)
	DAT <- data.frame(navec=rep(NA, point_num))
	names= c("NA_vec")
	for (lambda_var in lambda_vec){
		for (polynomial_m_var in polynomial_m_vec){
			filename = paste0('L = ', lambda_var, ' M= ', polynomial_m_var)
			message(filename)
			max_density <- max(retina_obj$azimuthal_data.datapoints[[1]]$z)
			tmp <- retinaplot( retina_obj, 
				return_fit=TRUE,
				spatial_res=spatial_res,
				contour_breaks_source  =  c(0,max_density), 
				col_breaks_source      =  c(0,max_density),
				col_levels=50,
				contour_levels=20,
				rotation=0,
				inner_eye_view=TRUE,
				lambda=lambda_var,
				polynomial_m=polynomial_m_var)
			fit_error_histogram(tmp[[1]], sub=filename, ylim=c(0,frequency_cap), xlim= error_x_limits)
			fit_object = tmp[[1]]
			DAT <- cbind(DAT, predictSE(fit_object))
			names <- c(names, filename)
		}
	}
	colnames(DAT) <- names
	return(DAT)
}

##' @title Merge Sampling Site location and Counting Frame counts.
##' @description
##' X represents the row index on the sampling grid, and Y represents the column index. Numbers start at 1.
##' @param location dataframe with a column for samplingsite, x and y. 
##' @param counts dataframe with a column for samplingsite and count.
##' @param ... further arguments passed to or from other methods.
##' @return xycount Dataframe of (x,y, count)
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family munge
##' @export
ssite_merge <- function(location, counts, ...){
	xycount = merge(location, counts, by = "samplingsite")
	xycount = subset(xycount, select=c('x','y','count'))
	return(data.frame(xycount))
}

##' @title ImageJ and Retistruction of Input Retinal Data from Flattened State
##' @description
##' Converts flat retina to a spherical representation of the data. Spherical coords path must contain:  outline.roi  Datapoints in barycentric coordinates. P.csv  Datapoints on reconstructed sphere in cartesian coordinates. T.csv  Datapoints on reconstructed sphere in spherical coordinates.  xyz.csv  Cell density site locations with density measurements. Columns of x(integer column starting at 1 at the left) , y (integer column starting at 1 on the bottom) and z (counts)
##' @param path Directory which contains the retistruct files.
##' @param height Height in microns of the stereologic counting frame
##' @param width Width in microns of the stereologic counting frame
##' @param IJ_limits data.frame with min, max and delta values for both X and Y with respect to the contour image in ImageJ coordiantes. For setting the calibration of the sampling sites to the sphere.
##' @param falciform (boolean) True by default, meaning there is a file called falc.txt within the path.
##' @return data.frame with phi(latitude), lambda(longitude) and Z (cells per square millimeter).
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @references https://r-forge.r-project.org/scm/?group_id=1436
##' @import retistruct
##' @export
spherical_coords <- function(path, height, width, IJ_limits, falciform=TRUE){
	#Read in the xyz dataset from stereology data collection.
	xyz <- import_xyz(path)
	#Read in the Falciform Process x y outline coordinates.
	if (falciform==TRUE){
		falc <- read.table(paste0(path, "/falc.txt"), col.names=c("x", "cyan")) 
	} else {
		#If not, just put in NULL.
		falc <- data.frame()
	}
  

	# Import Counting Frame surface area in microns^2
	xy_IJ_cols <- coordinate_IJ(xyz,  
						  IJ_limits$maxX,
						  IJ_limits$maxY,
						  IJ_limits$minX,
						  IJ_limits$minY,
						  IJ_limits$deltaX,
						  IJ_limits$deltaY,
						  inversion=FALSE)
	# Save output/species/datapoints.csv, create folder if doesnt exist #saves only the x,y coordinates
	DAT <- data.frame(xy_IJ_cols[,1], xy_IJ_cols[,2])
	xyz_len <- length(xy_IJ_cols[,1])
	# Append the falciform process to the end of the dataset.
	colnames(DAT) <- c("x","cyan")
	DAT<-rbind(DAT,falc) #Tag the falciform process onto the end
	combined_len <- length(DAT[,1])
	#Save the points as datapoints.csv. This includes XYZ datapoints and the falciform process.
	write.csv(DAT,file=paste0(path, "/datapoints.csv"), row.names=FALSE) 

	#Post-image_J_markup
	radian_data <- dss_retistruct_processing(path)
	dss_object<- getDss(radian_data)
	dss <- getDssRemoved(radian_data)
	#Extract the density measurement and falciform outline coordinates from the 'datapoints' dss set.
	falciform_outline<- dss$x[((xyz_len):combined_len),]
	falciform_outline<- falciform_outline[(2:length(falciform_outline[,1])),]
	density_locations<- dss$x[1:xyz_len,]
	dss_coords <- degrees(data.frame(density_locations))
	falciform_outline <- degrees(data.frame(falciform_outline))
	#Combine coordinates with sampling site densities
	#Get rid of the NA points http://stackoverflow.com/questions/4862178/remove-rows-with-nas-in-data-frame
	data_w_NA <- cbind(dss_coords,z = xyz[,3])
	trimmed_data <- data_w_NA[complete.cases(data_w_NA[,1:2]),] 
	falciform_outline <- falciform_outline[complete.cases(falciform_outline[,1:2]),] 
	#convert to density per square millimeter
	trimmed_data$z <- unlist(lapply(trimmed_data$z, 
		function(x) count_to_rho(x, height=height, width=width)))
	return(list(trimmed_data=trimmed_data, falciform_outline=falciform_outline))
}

##' @title ImageJ Coordinate Conversion
##' @param RET_count_data XYZ dataframe. X and Y are integer values starting at 1.
##' @param maxX The X value of the furthest sampling site at the far right of the picture.
##' @param maxY The Y value of the furthest sampling site at the top of the picture
##' @param minX The X value of the furthest sampling site at the far left of the picture
##' @param minY The Y value of the furthest sampling site at the bottom of the picture
##' @param deltaX Average change in pixels in the X axis between two sampling sites.
##' @param deltaY Average change in pixels in the Y axis between two sampling sites.
##' @param inversion Whether or not to invert
##' @param ... arguments passed to or from other methods.
##' @return RET_count_data Retinal count data in ImageJ coordinates. $x, $y and $z
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
coordinate_IJ <- function(RET_count_data,	maxX,
											maxY,
											minX,
											minY, 
											deltaX, 
											deltaY, inversion,...)
{
	#data = XY sampling grid integer values are converted into imageJ coordinates based on image.
	# make sure that maxY and minY are inputted as negative values.
	if (min(RET_count_data$x)==1){
		message("x starts at 1, no change committed")}
	if (min(RET_count_data$y)==1){
		message("y starts at 1, no change committed")}
	if (min(RET_count_data$x)!=1){
		offset<-min(RET_count_data$x)-1
		message(paste0(offset, " is the offset fixed for x"))
		RET_count_data$x <- RET_count_data$x-offset
		message("x doesn't start at 1, change committed")
	}
	if (min(RET_count_data$y)!=1){message("y doesn't start at 1, change committed")
		offset<-min(RET_count_data$y)-1
		message(paste0(offset, " is the offset fixed for y"))
		RET_count_data$y <- RET_count_data$y-offset
	}
	RET_count_data$x <- minX + (RET_count_data$x-1)*deltaX
	if (inversion==TRUE){
		RET_count_data$y <- abs(maxY + (RET_count_data$y-1)*deltaY) #artificially made positive to meet IJ reqs
		message("INVERSION!, max added")
	}
	if (inversion==FALSE) {
		RET_count_data$y <- abs(minY + (RET_count_data$y-1)*deltaY) #artificially made positive to meet IJ reqs
		message("No inversion, min added")
	}
	  return(RET_count_data)
}

##' @title Retistruct Wrapper
##' @description
##' Reads in the folder, reconstructs, and returns the retinal object.
##' @param path string path to the folder with retinal data.
##' @return rad Retinal data in radian units.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}
##' @references Sterratt et. al. 2013
dss_retistruct_processing <- function (path){ 
	rad <- retistruct.read.markup(retistruct.read.dataset(path))
	rad <- retistruct.reconstruct(rad)## Reconstruct (computation intensive)
	return(rad)
}

##' @title Conversion from Radians to Degrees
##' @param radian_coords A data.frame with two columns- first is phi, second is lambda
##' @return data.frame Phi and lambda in a data.frame.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
degrees <- function(radian_coords){
	phi <- radian_coords[,1]*(180/pi) #convert to degrees, phi is out latitude
	lambda <- radian_coords[,2]*(180/pi) #convert to degrees, lambda is our longitude, use this to rotate about the eye through-axis
	return(data.frame(phi=phi,lambda=lambda))
}

##' @title Convert count to Density
##' @description
##' Uses cell count to compute density within the counting frame. Units are set to microns
##' @param count Number of retinal cells counted within the counting frame.
##' @param height Counting frame width in microns
##' @param width Counting frame width in microns
##' @return rho Density in units of retinal ganglion cells per square millimeter
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}
##' @examples
##' count_to_rho(10,25,25)
##' #16,000 cells/mm^2
##' @export
count_to_rho <- function(count, height, width){
	rho = count/(height*width*(1e-6))
	return(rho)
}


##' @title Spherical Plot visualization
##' Uses retistruct to create lat/lon coordinates
##' @param trimmed_data Three-column dataset with datapoint locations in latitude/longitude degree format.
##' @return Spherical visualization
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @import rgl sphereplot
##' @export
sphere_visualize<- function(trimmed_data){
	rgl.sphgrid(radius=1,longtype="D",deggap=30,col.lat="transparent",col.long="black") #make a quick sphere
	rgl.sphpoints(trimmed_data[,2], trimmed_data[,1], shininess=50.0, radius=1,size=6, deg=TRUE, col="blue4")#quick check to make sure it is covering the bottom hemisphere
}

##' @title Single Tps Fit Error Plot
##' @param x the fit object
##' @param main title of the plot, NULL by default
##' @param ... further arguments passed to or from other methods.
##' @return Error histogram in base R
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @export
fit_error_histogram <- function(x, main=NULL, ...){
	hist(predictSE(x), col="black", xlab="Fit Error (RGC/sq.mm)", ...)
}


##' @title Retinal Krig Fit Plots
##' @param x the fit object
##' @param digits set to 4 arbitrarily
##' @param which set arbitrarily to 1:4
##' @param ... further arguments passed to or from other methods.
##' @return Histograms and Error scatterplots in base R
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @references Fields Package
##' @export
fit_plots <- function(x, digits = 4, which = 1:4, ...) {
	par(mfrow=c(2,2))
	out <- x
	#
	#   don't do plots 2:4 if a fixed lambda
	#
	if (x$fixed.model) {
		which <- 1
	}
	fitted.values <- predict(out)
	std.residuals <- (out$residuals * sqrt(out$weights))/out$shat.GCV
	if (any(which == 1)) {
		temp <- summary(out)
		# plot.new(xlim = c(min(fitted.values), max(fitted.values)), ylim = c(min(out$y), max(out$y)) ylab = "Y", xlab = " Predicted density", 
		#     bty = "n", ...)
		# points(fitted.values, out$y, col=black, alpha = 0.5)
		# abline(0, 1)
		plot(fitted.values, out$y, ylab = "Y", xlab = " Predicted density", 
			bty = "n", ...)
		abline(0, 1)
		hold <- par("usr")
		text(hold[1], 0.95*range(out$y)[2], paste(" R**2 = ", format(round(100 * 
			temp$covariance, 2)), "%", sep = ""), cex = 0.8, 
			adj = 0)
	}
	if (any(which == 2)) {
		plot(fitted.values, std.residuals, ylab = "(STD) Residuals", 
			xlab = " Predicted density", bty = "n", ...)
		yline(0)
		hold <- par("usr")
		text(hold[1], 0.95*range(std.residuals)[2], paste(" RMSE =", format(signif(sqrt(sum(out$residuals^2)/(temp$num.observation - 
			temp$enp)), digits))), cex = 0.8, adj = 0)
	}

	if (any(which == 4)) {
		hist(std.residuals, main="", col="black", breaks = 10, xlab = "STD Residual")
	}
	message("Predicting Fit Error")
	hist(predictSE(x), col="black", xlab="Fit Error (RGC/sq.mm)", main=NULL)
	return(predictSE(x))
}
### Set contour_breaks based on requested source
##' @title Set Contour Breaks Based on Requested Source
##' @description This function will make a set of contour topography lines. Code from http://stackoverflow.com/questions/10856882/r-interpolated-polar-contour-plot was highly modified to meet retinal plotting funtionality.
##' @param contour_breaks_source See fit_plot_azimuthal
##' @param z See fit_plot_azimuthal
##' @param contour_levels See fit_plot_azimuthal
##' @param Mat See fit_plot_azimuthal
##' @value 0 returns 0 if no issue
define_contour_breaks <- function(contour_breaks_source, z, contour_levels, Mat) {
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
  return(contour_breaks)
}

### Add contours to the retina plot
##' @title Print contour lines onto the retina plot
##' @description Makes a set of contours to the retinaplot. Modified code from http://stackoverflow.com/questions/10856882/r-interpolated-polar-contour-plot was highly modified to meet retinal plotting funtionality.
##' @param minitics See fit_plot_azimuthal
##' @param Mat See fit_plot_azimuthal
##' @param contour_breaks See fit_plot_azimuthal
##' @import grDevices
##' @value 0 returns 0 if no issue
add_contours <- function(minitics, Mat, contour_breaks, xy){
	require(grDevices)
	CL <- grDevices::contourLines(x = minitics, y = minitics, Mat, levels = contour_breaks)
	A <- lapply(CL, function(xy){
  			lines(xy$x, xy$y, col = gray(.2), lwd = .5)
		})
	return(0)
}

##' @title Initiate the Square Matrix plot to prepare for polar plotting
##' @description instantiates the square plotting area Modified code from http://stackoverflow.com/questions/10856882/r-interpolated-polar-contour-plot was highly modified to meet retinal plotting funtionality.
##' @param zlim See fit_plot_azimuthal
##' @param col See fit_plot_azimuthal
##' @param Mat See fit_plot_azimuthal
##' @param minitics See fit_plot_azimuthal
##' @value 0 returns 0 if no issue
init_square_mat_plot <- function(Mat, zlim, minitics, col){
	Mat[which(Mat<zlim[1])]=zlim[1]
  	Mat[which(Mat>zlim[2])]=zlim[2]
  	image(x = minitics, y = minitics, Mat , useRaster = TRUE, asp = 1, axes = FALSE, xlab = "", ylab = "", zlim = zlim, col = col)
  	return(0)
}

##' @title Compute fit error at original data points
##' @description compute the error at a set of desired points
##' @author Brian Cohn
##' @param x the x coordinates of the original points that were smoothed
##' @param y the y coordinates of the original points that were smoothed
##' @param thin_plate_spline_object output object from fields::Tps
##' @value error A data frame of the error at each of the original points
compute_thin_plate_spline_error <- function(x,y,thin_plate_spline_object) {
	return(data.frame(x=x,y=y, se=predictSE(thin_plate_spline_object)))
}


##' @title Draw circle radians about the center
##' @description Draw N radius lines (circles)
##' @author Brian Cohn
##' @param number_of_circles the number of radius lines to plot
plot_circle_radial_lines <- function(number_of_circles, color_hex = "#66666650"){
	circle <- function(x, y, rad = 1, nvert = 500){
	  rads <- seq(0,2*pi,length.out = nvert)
	  xcoords <- cos(rads) * rad + x
	  ycoords <- sin(rads) * rad + y
	  cbind(xcoords, ycoords)
	}
	for (i in number_of_circles){
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


##' @title Plot falciform process
##' @description smooths and plots the falciform process
##' @author Brian Cohn
##' @param falciform_x numeric vector of x coordinates
##' @param falciform_y numeric vector of y coordinates
plot_falciform_process <- function(falciform_x, falciform_y){
	fc_smoothed <- spline.poly(cbind(falciform_x,falciform_y), 50, k=10)
	polygon(fc_smoothed[,1], fc_smoothed[,2], col=rgb(0, 0, 0,0.5), lty="solid", border="gray42")
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
##' @param interp.type depreciated
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
  ### Plotting data (in cartesian) - will be converted to polar space.
  x, y, z, 
  ### Plot component flags
  contours=TRUE,   # Add contours to the plotted surface
  legend=TRUE,        # Plot a surface data legend?
  axes=TRUE,      # Plot axes?
  points=TRUE,        # Plot individual data points
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
  interp.type = 1, #1 = Thin plate spline , null if you are inputting your own predicted surface
  lambda=0.001, xyrelief=1,tmp_input = NULL, plot_suppress=FALSE,
  compute_error=FALSE,
  falciform_coords=NULL,
  falc2=NA, 
  polynomial_m=NULL,
  ...){ 
  
  
  #Used only when interp.type = 2

  
  # interpolate the data
  if (interp.type == 1){
  
	minitics <- seq(-outer.radius, outer.radius, length.out = spatial_res)
	grid.list = list(x=minitics,y=minitics) #choose locations to predict at
	t <- Tps(cbind(x,y),z,lambda=lambda, m = polynomial_m) #computationally intensive
	tmp <- predictSurface(t,grid.list,extrap=extrapolate)
	Mat = tmp$z
  }
  else {stop("interp.type value not valid")}
  if (compute_error==TRUE){
  	error <- compute_thin_plate_spline_error(x,y,thin_plate_spline_object)
  } else {
  	error <- NULL
  }




  if (plot_suppress == TRUE){
		error<-NULL
		return(list(t,tmp, error))
  }


  #turn values outside of the circle to NA
  markNA <- matrix(minitics, ncol = spatial_res, nrow = spatial_res) 
  Mat[!sqrt(markNA ^ 2 + t(markNA) ^ 2) < outer.radius] <- NA 
  
  zlim <- define_color_breaks_based_on_source(col_breaks_source,z, Mat)

  init_square_mat_plot(Mat, zlim, minitics, col)
  
  if (contours){ add_contours(minitics, Mat,
  	contour_breaks=define_contour_breaks(contour_breaks_source, z, contour_levels, Mat), xy)}
  
  #Produce spline polygon for falciform 1
	plot_falciform_process(falciform_coords$x, falciform_coords$y)

  #Plug in the secondary plot if it is available
  if (!is.na(falc2)) plot_falciform_process(falc2$x, falc2$y)


  # add interpolated point if desired
  if (points){
	points(x+0.005,y,pch=4, cex=0.5, col="gainsboro")
	points(x,y,pch=4, cex=0.5, col="black")
   
  }

  # add radial axes if desired
  if (axes){ 
	# internals for axis markup
	RMat <- function(radians){
	  matrix(c(cos(radians), sin(radians), -sin(radians), cos(radians)), ncol = 2)
	}    
	
	# draw latitude markings
	if (missing(circle.rads)){
	  circle.rads <- pretty(c(0,outer.radius-.4))
	}
	plot_circle_radial_lines(circle.rads)

	# put on radial spoke axes:
	axis.rads <- c(0, pi / 6, pi / 3, pi / 2, 2 * pi / 3, 5 * pi / 6)
	r.labs <- c(90, 60, 30, 0, 330, 300)
	l.labs <- c(270, 240, 210, 180, 150, 120)
	
	for (i in 1:length(axis.rads)){ 
	  endpoints <- zapsmall(c(RMat(axis.rads[i]) %*% matrix(c(1, 0, -1, 0) * outer.radius,ncol = 2)))
	  segments(endpoints[1], endpoints[2], endpoints[3], endpoints[4], col = "#66666650")
	  endpoints <- c(RMat(axis.rads[i]) %*% matrix(c(1.1, 0, -1.1, 0) * outer.radius, ncol = 2))
	  lab1 <- bquote(.(r.labs[i]) * degree)
	  lab2 <- bquote(.(l.labs[i]) * degree)
	  text(endpoints[1], endpoints[2], lab1, xpd = TRUE, family = "Palatino")
	  text(endpoints[3], endpoints[4], lab2, xpd = TRUE, family = "Palatino")
	}

	axis(2, pos = -1.25 * outer.radius, at = sort(union(circle.rads,-circle.rads)), labels = NA)
	text( -1.26 * outer.radius, sort(union(circle.rads, -circle.rads)),sort(union(circle.rads, -circle.rads)), xpd = TRUE, pos = 2,family = "Palatino")
  }
  
#add legend
  if (legend){
		par(mai = c(1,1,1.5,1.5))
		image.plot(legend.only=TRUE, col=col, zlim=zlim, family = "Palatino")

  }

return(list(t,tmp, error))
}

##' @title Constructor for RecontructedDataset object *EDITED
##' @description This function was edited by Brian Cohn on 05/30/2014 in order to return DSS. One line was added.
##' @param r Object that of clases \code{reconstructedOutline} and
##' \code{dataset}.
##' @param report Function used to report progress.
##' @return \code{\link{ReconstructedDataset}} object containing the input
##' information and the following modified and extra information:
##' \item{\code{Dsb}}{Datapoints in barycentric coordinates}
##' \item{\code{Dsc}}{Datapoints on reconstructed sphere in cartesian coordinates}
##' \item{\code{Dss}}{Datapoints on reconstructed sphere in spherical coordinates}
##' \item{\code{Ssb}}{Landmarks in barycentric coordinates}
##' \item{\code{Ssc}}{Landmarks on reconstructed sphere in cartesian coordinates}
##' \item{\code{Sss}}{Landmarks on reconstructed sphere in spherical coordinates}
##' @author David Sterratt
##' @import geometry
##' @export
getDssRemoved <- function(r, report=message) {
  ###Function by David Sterratt, 2013
  report("Inferring coordinates of datapoints")
	#internal functions pasted in for convenience
  Dsb <- list() # Datapoints in barycentric coordinates
  Dsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Dss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ds) & (length(r$Ds) > 0)) {
	for (name in names(r$Ds)) {
	  Dsb[[name]] <- tsearchn(r$P, r$T, r$Ds[[name]])
	  oo <- is.na(Dsb[[name]]$idx)     # Points outwith outline
	  if (any(oo)) {
		warning(paste(sum(oo), name, "datapoints outwith the outline will be ignored."))
	  }
	  Dsb[[name]]$p   <- Dsb[[name]]$p[,,drop=FALSE]
	  Dsb[[name]]$idx <- Dsb[[name]]$idx
	  # Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]]) #ORIGINAL LINES
	  # Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
	  Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]]) #MODIFIED LINE, BRIAN COHN 05/30/2014
	  Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
	}
  }
  return(Dss)#ADDED LINE, BRIAN COHN 05/30/2014
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



##' @title Counterclockwise rotation about the origin
##'
##' @param x vector of x coordinates (int/float)
##' @param y vector of y coordinates (int/float)
##' @param theta (float|int) angle in degrees.
##' @return xy Data.Frame with $x and $y rotated coordinates
##' 
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family trig_fns
##' @examples
##' x<- c(1,2,3,4,5)
##' y<- c(3,3,3,3,3)
##' theta <- 30
##' cartesian_rotation(x,y,theta)
##' @export
cartesian_rotation<- function (x, y, theta){
  deg <- function(radians) 180*radians/pi
  rad <- function(degrees) degrees*pi/180
  theta = rad(theta)
  x_new <- x*cos(theta) - y*sin(theta)
  y_new <- y*cos(theta) + x*sin(theta)
  return(data.frame(x=x_new, y=y_new))
}
 

##' @title Add Degrees
##' @param theta Angle in degrees
##' @param degrees Number of degrees to add to theta
##' @return newtheta Adjusted between -180 and 180 degrees.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family trig_fns
##' @examples
##' add_degrees(30, 180) #you should get -150
##' @export
add_degrees <- function(theta, degrees){
	if (theta+degrees > 180){
		theta <- (theta+degrees - 360)
	} else if (theta+degrees < -180) {
		theta <- (theta+degrees + 360)
	} else {
		theta <- (theta+degrees)
	}
	return(theta)
}




##' @title Rotation error measurement for two similar retinal maps.
##' @description
##' Rotates one map 360 degrees, and for every value of rotation between -180 to 180, the 
##' mean absolute difference between the maps' rho values is recorded
##' @param map1 retinal data: This map will not be rotated, and will serve as the reference. 
##' @param map2 retinal data: This map will be rotated and superimposed upon map1.
##' Has components $x, $y and $z in azimuthal projection coordinates.
##' @param ... further arguments passed to or from other methods.
##' reference map. latitude in [,1], longitude in [,2], rho in [,3].
##' @param spatial_res Defines the number of pixels across the interpolation grid for every comparison. Default is 16 for speed.
##' @param theta_interval Defines the interval at which the rotation value will be changed when moving to the next comparison. Default is 10 in units of degrees.
##' @return df_spin Matrix of errors (in second column) for a given map rotation.
##' @import mapproj
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family species_average
##' @export
rotation_optimize <- function(map1, map2, spatial_res = 32, theta_interval = 30, projection_type="azequidistant",...){
  #define prediction frame for both maps
  map1 <- map1$azimuthal_data.datapoints[[1]]
  map2 <- map2$trimmed_data
  minitics <- seq(-1.6, 1.6, length.out = spatial_res)
	grid.list = list(x=minitics,y=minitics)

	fit1<- Tps(cbind(map1$x,map1$y),map1$z,...)
	map1surface <- predictSurface(fit1, grid.list, extrap=FALSE)$z
  #initialize df_spin dataframe to hold the RMSE values of each rotation value
  
	theta <- seq(-180, 180, by=theta_interval)
  df_spin <- matrix(nrow=length(theta),ncol=2)

  #loop through spin theta
  for (rotation in 1:length(theta)){
	#Create an azimuthal equidistant map projection with "rotation" variable.
	  az<-mapproject(x=map2[,2], map2[,1], 
		 projection=projection_type, orientation=c(-90,0,theta[rotation])) 
	  az$z <- map2[,3]
	#get prediction frame
	suppressWarnings(fit2 <- Tps(cbind(az$x,az$y),az$z,...))
	map2surface <- predictSurface(fit2, grid.list, extrap=FALSE)$z
	map_difference <- map2surface - map1surface #matrix difference
	  meandiff <- mean(abs(map_difference), na.rm=TRUE)
	  df_spin[rotation,] <-  c(theta[rotation], meandiff)
  }
  # optimal_rotation(df_spin)
  colnames(df_spin) <- c("degree", "abs_mean_diff")
  df_spin[,1] <- unlist(lapply(df_spin[,1], add_degrees, degrees=90))
  return(df_spin)
}



##' Optimal rotation (for min_err) for rotation error matrix
##' Computes the argmax of the min error, returning both the error and its theta rotation value in degrees.
##' @param df_spin Two column matrix of theta (rotation in degrees) at [,1] and error at [,2].
##' @param quiet logical. when false, it will print what the absolute mean difference is between maps.
##' @return optimal_degree The optimal (minimum) degree and its error value.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family species_average
##' @export
optimal_rotation<- function(df_spin, quiet=TRUE){
  optimal_degree <- df_spin[which.min(df_spin[,2]), ] 
  if (!quiet) {
  	message(paste("Absolute mean difference between maps (", optimal_degree[2], "RGC/sqmm)  minimized when map2 ccw is ", optimal_degree[1], "degrees"))
  }
  return(optimal_degree)

}

##' Reflect matrix across vertical axis
##' @param mat a matrix (doesn't have to be square)
##' @return matrix reflected matrix
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu} Lars Schmitz
##' @family internal
reflect_across_vertical_line <- function(mat){
	num_cols <- length(mat[1,])
	num_rows <- length(mat[,1])
	for (i in 1:num_rows) {
		mat[,i] <- rev(mat[,i])
	}
	return(mat)
}

##' @title Make a Composite Map
##' @param map1 Fixed retina object
##' @param map2 retina object subject to rotation.
##' @param rotation Boolean, whether or not spin optimization is used.
##' @param plot_rotation Boolean, whether or not spin optimization plot is shown.
##' @param spatial_res Default is 16. A (spatial_res*spatial_res) matrix will be created for spin optimization comparisons and compose the average.
##' @param ... further arguments passed to or from other methods.
##' @return plot Composite Plot data
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family species_average
##' @export
composite_map <- function (map1, 
						   map2, 
						   rotation=TRUE, spatial_res, plot_rotation=FALSE, projection_type="azequidistant",
						   ...){

	###STEP 1 Convert map2 data into an azimuthal equidistant plot projection
		#map1
		x1 <- map1$azimuthal_data.datapoints[[1]]$x
		y1 <- map1$azimuthal_data.datapoints[[1]]$y
		z1 <- map1$azimuthal_data.datapoints[[1]]$z
	#map2
		x2 <- map2$azimuthal_data.datapoints[[1]]$x
		y2 <- map2$azimuthal_data.datapoints[[1]]$y
		z2 <- map2$azimuthal_data.datapoints[[1]]$z

	# par(mfrow=c(2,1))
	# retinaplot(map1, spatial_res=50)
	# retinaplot(map2, spatial_res=50)
	# par(mfrow=c(1,1))

	theta <- -90 #default value is -90 (no rotation)
	if (rotation){
		message("Optimizing Rotation")
		#identify error at each rotation
		rotation_df <- rotation_optimize(map1, map2, 64)
		if (plot_rotation) {
		 	plot_rotation_optimize(rotation_df)
		}
		#find optimal rotation for map2 with respect to map1
		theta <- optimal_rotation(rotation_df)
		ori <- c(-90, 90, theta[1])
			az <- mapproject(x=map2$trimmed_data[,2], map2$trimmed_data[,1], 
						   projection=projection_type, orientation=ori)
		x2 <- az$x
		y2 <- az$y
	}
	message("Fitting Tps Model for Map1")
	##Step 2 Interpolate map1 and map2 azimuthal data with the same interpolation options
	map1fit <- fit_plot_azimuthal(x1,y1,z1, spatial_res,
									plot_suppress=TRUE, extrapolate=TRUE,
									outer.radius=pi/2,
									falciform_coords = map1$azimuthal_data.falciform[[1]],...)
	message("Fitting Tps Model for Map2")
	map2fit <- fit_plot_azimuthal(x2,y2,z2, spatial_res,
									plot_suppress=TRUE, extrapolate=TRUE,
									outer.radius=pi/2,
									falciform_coords = map2$azimuthal_data.falciform[[1]],...)
	maps <- list(map1fit, map2fit)
	MAT1 <- maps[[1]][[2]]$z	
	MAT2 <- maps[[2]][[2]]$z	
	#90 degree CCW rotation
	# reflect across a vertical line in the center of the retina. This changes the view to an inner view
	MAT1 <- reflect_across_vertical_line(MAT1)
	MAT2 <- reflect_across_vertical_line(MAT2)

	# MAT1 <- sapply(nrow(MAT1):1, function(i) MAT1[i, ])
	# MAT2 <- sapply(nrow(MAT2):1, function(i) MAT1[i, ])
	# map_composites(MAT1,MAT2, ...)
	composite_pred <- (MAT1 + MAT2)/2.0
	return (composite_pred)
}


##' @title Composite two Retinal Map Matrices
##' @param MAT1 Matrix of RGC densities
##' @param MAT2 Matrix of RGC densities
##' @param col_levels Number of levels to plot (if showing plot.)
##' @param showplots False by default.
##' @param ... further arguments passed to or from other methods.
##' @return MAT_composite Square matrix: Composite Map (adjusted to max,min of each map, with new peak equal to the mean of map1 and map2's peaks)
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu}, Lars Schmitz
##' @family species_averags
##' @export
map_composites <- function(MAT1,MAT2, col_levels=5, showplots=FALSE,...){
	#Matrix Combinations
	max1 <- max(MAT1, na.rm=TRUE)
	max2 <- max(MAT2, na.rm=TRUE)
	avMAX <- mean(max1,max2)

	min1 <- min(MAT1, na.rm=TRUE)
	min2 <- min(MAT2, na.rm=TRUE)
	avmin <- mean(min1,min2)
	range1 <- (max1-min1)
	range2 <- (max2-min2)
	SUM <- MAT1+MAT2
	DIFF <- MAT1-MAT2
	MEAN <- (MAT1+MAT2)/2.0
	MAT1_rel <- (MAT1-min1)/range1
	MAT2_rel <- (MAT2-min2)/range2
	rel_mean <- (MAT1_rel+MAT2_rel)/2.0
	av_adjusted <- rel_mean*avMAX
	message("plotting")
	if (showplots == TRUE){
		plot.new()
		par(mfrow=c(3,3))
		#Simple
		image.plot(MAT1, col=heat.colors(col_levels), main="map1", useRaster=TRUE)
		image.plot(MAT2, col=heat.colors(col_levels), main="map2", useRaster=TRUE)
		#Combinations
		image.plot(SUM, col=heat.colors(col_levels), main="Sum", useRaster=TRUE)
		image.plot(DIFF, col=heat.colors(col_levels), main="Difference", useRaster=TRUE)
		image.plot(MEAN, col=heat.colors(col_levels), main="Mean", useRaster=TRUE)
		#Relative Maps
		image.plot(MAT1_rel, col=heat.colors(col_levels), main="map1 relative", useRaster=TRUE)
		image.plot(MAT2_rel, col=heat.colors(col_levels), main="map2 relative", useRaster=TRUE)
		image.plot(rel_mean, col=heat.colors(col_levels), main="Relative density mean", useRaster=TRUE)
		image.plot(av_adjusted, col=heat.colors(col_levels), main="Av, adj. to mean peak RGC density", useRaster=TRUE)
	}
	return(av_adjusted)
}

##' @title Plot rotation optimization scatter plot with minima.
##' @description
##' A base-graphics plot (rotation along the x axis,
##' Absolute Mean difference on the y axis) of the rotation optimization values, and a vertical line denoting the optimized rotation value.
##' @param rotation_df data.frame from spin optimization
##' @export
plot_rotation_optimize<- function (rotation_df){
	plot(rotation_df, xlim= c(-181,181), pch=20, 
		ylim=c(min(rotation_df[,2]) - 1000, max(rotation_df[,2])),
		xlab= "Rotation (ccw) [Degrees]", 
		main="Mean Absolute Difference",
		ylab=expression(mu(group("|", "map1-map2", "|"))))
	abline(v=optimal_rotation(rotation_df), lty="dotted")
}



#Utility functions:

##' @title Distance between the range
##' @param x vector of numbers
##' @return the difference between the limits of the range.
range_len <- function(x){
	limits <- range(x)
	diff <- limits[2]-limits[1]
	return (diff)
}
##' @title Reorder columns by function
##' @param X data.frame
##' @param FN A function that takes in a vector and returns one value.
##' @return X_new data.frame of columns ordered in decreasing value from left to right according to each function(column)
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu} Lars Schmitz
reorder_columns <- function(X, FN){
	Xcol_fn_applied <- apply(X, 2, FN)
	X_new <- X[,c(names(sort(Xcol_fn_applied, decreasing=TRUE)))]
	return(X_new)
}
