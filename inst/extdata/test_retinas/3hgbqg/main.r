main_3hgbqg <- function(path_to_diagram_retina){

  require(retina)

  #parameters
  LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
  RESOLUTION_var <- 500 #Plot width in pixels


  #defining coordinates
  IJ <-data.frame(maxX = 1429,
                  maxY = -85,
                  minX = 108,
                  minY = -989,
                  deltaX = (1429-108)/21,
                  deltaY = (989-85)/14)


  #assembling the retina object
  my_retina <- retina_object(
      path = path_to_diagram_retina,

      #Eye Measurements from dissection
          ED = 5,    #Eye diameter (mm)
          AL = 5,   #Eye axial length (mm)
          LD = 2,    #Eye lens diameter (mm)

      #Stereology Parameters
          height = 30 ,# height of the counting frame in microns
          width  = 30, # width of the counting frame in microns

      #Plotting Parameters
          lambda = LAMBDA_var,          #see fields::Tps for more information
          extrapolate = TRUE ,          #Predicts densities to the equator.
          spatial_res = RESOLUTION_var,
          rotation_ccw = -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
          plot_suppress=TRUE,
      #ImageJ Datapoint Calibration Measurements
          IJcoords = IJ
          )


  #plot

  retinaplot(my_retina)
}
