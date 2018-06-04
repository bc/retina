main_8pm223_direct_csv_coords <- function(path_to_diagram_retina){

LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels

my_retina <- retina_object(
  path = path_to_diagram_retina,

  #Eye Measurements from dissection
  ED = 4.8,    #Eye diameter (mm)
  AL = 3.43,   #Eye axial length (mm)
  LD = 1.8,    #Eye lens diameter (mm)

  #Stereology Parameters
  height = 25 ,# height of the counting frame in microns
  width  = 25, # width of the counting frame in microns

  #Plotting Parameters
  lambda = LAMBDA_var,          #see fields::Tps for more information
  extrapolate = TRUE ,          #Predicts densities to the equator.
  spatial_res = RESOLUTION_var,
  rotation_ccw = -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
  plot_suppress=TRUE)

retinaplot(my_retina)
return(my_retina)
}
