main_gv3igs <- function(path_to_diagram_retina){
#parameters
LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels


#defining coordinates
IJ <-data.frame(maxX = 1084,
                maxY = -36,
                minX = 100,
                minY = -1022,
                deltaX = (1084-100)/21,
                deltaY = (1022-36)/14)


#assembling the retina object
my_retina <- retina_object(
    path = path_to_diagram_retina,

    #Eye Measurements from dissection
        ED = 4.8,    #Eye diameter (mm)
        AL = 3.43,   #Eye axial length (mm)
        LD = 1.8,    #Eye lens diameter (mm)

    #Stereology Parameters
        height = 35 ,# height of the counting frame in microns
        width  = 35, # width of the counting frame in microns

    #Plotting Parameters
        lambda = LAMBDA_var,          #see fields::Tps for more information
        extrapolate = TRUE ,          #Predicts densities to the equator.
        spatial_res = RESOLUTION_var,
        rotation_ccw = -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
        plot_suppress=TRUE,
    #ImageJ Datapoint Calibration Measurements
        IJcoords = IJ)


#plot

# retinaplot(my_retina,rotation=270)
return(my_retina)
}
