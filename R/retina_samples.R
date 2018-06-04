##' 3hgbqg
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
main_3hgbqg <- function(path_to_diagram_retina){

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
  return(my_retina)
}


##' main_3hgbqg_direct_csv_coords
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
main_3hgbqg_direct_csv_coords <- function(path_to_diagram_retina){

  #parameters
  LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
  RESOLUTION_var <- 500 #Plot width in pixels

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
          plot_suppress=TRUE
          )

  #plot
  return(my_retina)
}


##' main_40oik5
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
main_40oik5 <- function(path_to_diagram_retina){
LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels

IJ <-data.frame(maxX = 944,
                maxY = -32,
                minX = 41,
                minY = -992,
                deltaX = 56,
                deltaY = 56)

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
  rotation_ccw = 0, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
  plot_suppress=TRUE,
  #ImageJ Datapoint Calibration Measurements
  IJcoords = IJ)


return(my_retina)
}


##' main_8pm223
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
main_8pm223 <- function(path_to_diagram_retina){

LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels

IJ <-data.frame(maxX = 981,
                maxY = -38,
                minX = 31,
                minY = -989,
                deltaX = 56,
                deltaY = 56)

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
  plot_suppress=TRUE,
  #ImageJ Datapoint Calibration Measurements
  IJcoords = IJ)


return(my_retina)
}




##' main_8pm223_direct_csv_coords
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
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


return(my_retina)
}

##' main_gv3igs
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
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


##' main_bjii4t
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
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


##' main_jhvbwt
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
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


return(my_retina)
}


##' main_raxp91
##' sample retina
##' @param path_to_diagram_retina path to folder
##' @return r result from retina_object
##' @export
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

return(my_retina)
}
