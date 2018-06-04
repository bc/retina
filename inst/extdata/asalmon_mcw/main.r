#parameters
LAMBDA_var <- 0.1 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #Plot width in pixels

#TODO change path. make sure this main.r script is just outside the diagram_retina folder.
path_to_retina_data_folder <- "~/Downloads/asalmon_mcw/diagram_retina"
outline_coordinates <- tear_markup_plot(path_to_retina_data_folder)
#todo verify these. see the tutorial here: https://github.com/bc/retina/blob/master/tutorial.md search for the part that starts with "#[V point of the tear,"
tear_coordinates_dataframe <- rbind(
    c(59,52,67),
    c(79,73,85),
    c(34,31,37))
                                                                    )
assemble_tear_file(tear_coordinates_dataframe, path_to_retina_data_folder)
assemble_point_coordinates_file(outline_coordinates, path_to_retina_data_folder)
# TODO fix nasal outline index, or set dorsal, set left or right eye
assemble_markup_file('left', path_to_retina_data_folder, nasal_outline_index=1, dorsal_outline_index=NA)


#assembling the retina object
my_retina <- retina_object(
    path = path_to_retina_data_folder,

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
        plot_suppress=TRUE)

retinaplot(my_retina, lambda=0.001, polynomial_m=   2, rotation=0)\
par(mfrow=c(2,2))
fit_plots(my_retina$fit_data1)
