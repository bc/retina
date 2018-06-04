main_bjii4t <- function(path_to_diagram_retina){

LAMBDA_var <- 0.1 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels


#defining coordinates
IJ <-data.frame(maxX = 1210,
                maxY = -46,
                minX = 72,
                minY = -928,
                deltaX = (1210-72)/20,
                deltaY = (928-46)/15)


#assembling the retina object
my_retina <- retina_object(
    path = path_to_diagram_retina,

    #Eye Measurements from dissection
        ED = 8.5,    #Eye diameter (mm)
        AL = 8,   #Eye axial length (mm)
        LD = 5,    #Eye lens diameter (mm)

    #Stereology Parameters
        height = 45 ,# height of the counting frame in microns
        width  = 45, # width of the counting frame in microns

    #Plotting Parameters
        lambda = LAMBDA_var,          #see fields::Tps for more information
        extrapolate = TRUE ,          #Predicts densities to the equator.
        spatial_res = RESOLUTION_var,
        rotation_ccw = -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
        plot_suppress=TRUE,
    #ImageJ Datapoint Calibration Measurements
        IJcoords = IJ)


#plot

retinaplot(my_retina)

findPeak <- function(my_retina){
	coord <- my_retina$trimmed_data
	peak_temp <- coord[which.max(coord$z),]
	peak <- data.frame(latitude=abs(peak_temp[,1]), longitude=peak_temp[,2], peak=peak_temp[,3])
	return(peak)
	}

message("PEAK DETERMINATION")
message(findPeak(my_retina))
return(my_retina)
#optional arguments
  #lambda=0.001
  #polynomial_m=2
  #rotation=30

#e.g. rotate map by 30 degrees: retinaplot(my_retina, rotation=30)
}
