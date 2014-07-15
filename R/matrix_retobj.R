##' @title Extract a matrix from a retinal object
##' @description
##' You can decide how much to rotate the retina before returning the matrix. You can also pass arguments to the polar plotter.
##' @param retina_object A retina
##' @param rotation in degrees. By default 0; set to zero rotation.
##' @param n the number by which to divide the matrix by. Useful when generating partial averages.
##' @param reflect logical, by default it is TRUE, in order to view the retina from outside the eye (rather than from within.)
##' @param ... Other arguments passed from other functions
##' @return mat A square matrix of the evaluated retina fit.
##' @author Brian Cohn \email{brian_cohn14@@pitzer.edu} and Lars Schmitz
##' @export
mat_from_ret_obj<- function(retina_object, rotation=0, n=1, reflect=TRUE, ...){
	x2 <- retina_object$azimuthal_data.datapoints[[1]]$x
	y2 <- retina_object$azimuthal_data.datapoints[[1]]$y
	z2 <- retina_object$azimuthal_data.datapoints[[1]]$z

	# Rotation
	ori <- c(-90, 90, rotation)
	az <- mapproject(x=retina_object$trimmed_data[,2], retina_object$trimmed_data[,1], 
					   projection="azequidistant", orientation=ori)
	x2 <- az$x #overwirte old value
	ifelse(reflect,x2<-(-1.0*x2), message(''))
	y2 <- az$y #overwrite old value

	# Fitting models
	map2fit <- fit_plot_azimuthal(x2,y2,z2,
									plot_suppress=TRUE, extrapolate=TRUE,
									outer.radius=pi/2.0,
									falciform_coords = retina_object$azimuthal_data.falciform[[1]],...)
	MAT <- map2fit[[2]]$z	
	#reflect the eye
	return(MAT/n)
}