main_jhvbwt <- function(path_to_diagram_retina){
#parameters
LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels


#defining coordinates
IJ <-data.frame(maxX = 1189,
                maxY = -81,
                minX = 26,
                minY = -986,
                deltaX = (1189-26)/19,
                deltaY = (986-81)/15)


#assembling the retina object
my_retina <- retina_object(
    path = path_to_diagram_retina,

    #Eye Measurements from dissection
        ED = 7,    #Eye diameter (mm)
        AL = 6.44,   #Eye axial length (mm)
        LD = 3.78,    #Eye lens diameter (mm)

    #Stereology Parameters
        height = 25 ,# height of the counting frame in microns
        width  = 25, # width of the counting frame in microns

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
}
