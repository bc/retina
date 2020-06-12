source("helper_functions.r")
temp_dir <- tempdir(check = TRUE)

roi_path <- "sample_retina/outline.roi"

# Coordinates of datapoints collected
measurements <- fread("sample_retina/datapoints.csv")
falciform_coords_raw <- fread("sample_retina/falciform.csv")
# Add values at each datapoint


#use this plot to define the tears
save_outline_indices_plot(load_roi(roi_path),measurements,"output/outline_coordinates.pdf")

# A tear is a row: c(middle tear point, before, after). Add commas as necessary
tear_coordinates_dataframe <- rbind(
	c(34,31,35),
	c(5,4,6),
	c(14,11,16),
	c(25,23,27)
	)
markup_information <- data.frame(eye_side=as.character("right"), nasal_outline_index=NA, dorsal_outline_index=8,phi0=0)
markup_information$eye_side <- as.character(markup_information$eye_side)

assemble_markup_file(markup_information$eye_side, temp_dir, nasal_outline_index=markup_information$nasal_outline_index, dorsal_outline_index=markup_information$dorsal_outline_index, phi0=markup_information$phi0)
# Create a temporary directory to hold onto the incremental work
file.copy(roi_path, file.path(temp_dir,"outline.roi"),overwrite=TRUE)
tear_df <- assemble_tear_file(tear_coordinates_dataframe, temp_dir)

r_outline <- retistruct.read.dataset(temp_dir) # Expect scale bar warning
r_markup <- retistruct.read.markup(r_outline) # make sure there is no P.csv already in the folder
r_reconstructed <- retistruct.reconstruct(r_markup)
# projection(r_reconstructed)

grid_coordinates <- grid_within_bounding_box(load_roi(roi_path)$coords,100)
grid_dt <- generate_projection_data(grid_coordinates, r_reconstructed)
landmarks <- directional_landmarks(grid_dt)

outline_coordinates <- extract_wholemount_outline(load_roi(roi_path))
outline_dt <- generate_projection_data(outline_coordinates, r_reconstructed)

measurement_coordinates <- apply_retistruct_inversion_to_datapoints(measurements,load_roi(roi_path))
measurement_dt <- generate_projection_data(measurement_coordinates, r_reconstructed)

falciform_coordinates <- apply_retistruct_inversion_to_datapoints(falciform_coords_raw,load_roi(roi_path))
falciform_dt <- generate_projection_data(falciform_coordinates, r_reconstructed)


## Flatmount plotting

p <- ggplot() + coord_fixed() + theme_classic()
p <- p + geom_polygon(aes(x,y), data = outline_dt, col="black", alpha=0) # show outline
p <- p + geom_point(aes(x, cyan, col=measurement), data = measurement_dt)  + scale_colour_gradient2()
p <- p + geom_label(aes(x,y,label=name),data=landmarks, alpha=0.5)
p <- p + geom_polygon(aes(x,cyan),data=falciform_dt, alpha=0, col="black")
ggsave("output/flatplot.pdf",p,width=8, height=8)


## Reconstructed plotting of the hemisphere

#remove the NA's where there were no points
dtt_projected <- na.omit(measurement_dt, cols="azi_x")

# PDF output is cleaner
# pdf("output/my_retina_figure.pdf", width=11,height=8.5, useDingbats=FALSE)
png("output/my_retina_figure.png", width=1200,height=800)
fit_data <- fit_plot_azimuthal(dtt_projected$azi_x,
	dtt_projected$azi_y,
	z = dtt_projected$measurement,
	outer_radius = 1.6,
	spatial_res = 1000,

    lambda = 0.01,
	col_levels = 50,
	contour_levels = 20,
	extrapolate = FALSE,

    compute_error = TRUE,
	eye_diameter = 0,
	axial_len = 0,
	falciform_coords = NA)
polygon(falciform_dt$azi_x,falciform_dt$azi_y, col="black")
dev.off()



