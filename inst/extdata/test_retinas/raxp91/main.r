main_raxp91 <- function(path_to_diagram_retina){
#parameters
LAMBDA_var <- 0.1 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels


#defining coordinates
IJ <-data.frame(maxX = 1016,
                maxY = -48,
                minX = 25,
                minY = -923,
                deltaX = (1016-25)/18,
                deltaY = (923-48)/17)


#assembling the retina object
my_retina <- retina_object(
    path = path_to_diagram_retina,

    #Eye Measurements from dissection
        ED = 8.5,    #Eye diameter (mm)
        AL = 7,   #Eye axial length (mm)
        LD = 4.5,    #Eye lens diameter (mm)

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

retinaplot(my_retina, lambda=0.001, polynomial_m=	2, rotation=30)
return(my_retina)
}
