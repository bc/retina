source('process_retinas.R')
# If spin_optimization is disabled, no rotation will be performed.
# composite_map evaluates the density across an interpolation grid at the desired *spatial_res = n* spatial resolution, and saves this data into an n by n matrtix.
# Finally, both matrices are added, and divided by two to produce the 'matrix average'. This map can be plotted by passing the matrix to the plot_from_MAT function.
# **note: the resultant map does not contain a combined falciform process, and the original sampling site points are not visible on the species average.
RESOLUTION_var <- 512
Pmol_av <- composite_map(Pmol_752, Pmol_753, rotation=TRUE, spatial_res=RESOLUTION_var)

# Pmol_av <- composite_map(Pmol_753, Pmol_752, rotation=TRUE, spatial_res=RESOLUTION_var)
boxplot(Pmol_752$azimuthal_data.datapoints[[1]]$z,
		Pmol_av,
		Pmol_753$azimuthal_data.datapoints[[1]]$z ,
		names=c("Sample 752","Average","Sample 753"),
		xlab="Retinal ganglion cells per square mm",
		horizontal=TRUE,
		pch=20, cex=0.5, col="lightgray")
rho_max <- max(Pmol_av)
rho_min <- min(Pmol_av)
#plot the species composite for Pmol
retinaplot(Pmol_752,
	contour_breaks_source  =  c(rho_min,rho_max),
	col_breaks_source      =  c(rho_min,rho_max),
		col_levels=50,
	contour_levels=20,
	spatial_res=RESOLUTION_var)
retinaplot(Pmol_753,
	contour_breaks_source  =  c(rho_min,rho_max),
	col_breaks_source      =  c(rho_min,rho_max),
		col_levels=50,
	contour_levels=20,
	spatial_res=RESOLUTION_var)
retinaplot(Pmol_753,
	contour_breaks_source  =  c(rho_min,rho_max),
	col_breaks_source      =  c(rho_min,rho_max),
		col_levels=50,
	contour_levels=20,
	spatial_res=RESOLUTION_var, rotation=-30)
plot_from_MAT(Pmol_av,
	extrapolate=TRUE,
	spatial_res = RESOLUTION_var,
	col_levels=50,
	contour_levels=20,
	contour_breaks_source  =  c(rho_min,rho_max),
	col_breaks_source      =  c(rho_min,rho_max),
	z1 = Pmol_753$azimuthal_data.datapoints[[1]]$z,
	z2 = Pmol_752$azimuthal_data.datapoints[[1]]$z
	)
