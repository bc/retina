par(mfrow=c(1,1))
# pdf("analysis_demo.pdf", height=7, width=8.1, useDingbats=FALSE, family="Palatino")

#Lets make a species average.
# Use the *composite_map* function to make a species average of two retinae.
# The user can enable "spin optimization", which we designed to help users overlay maps with minimal orientation error.
# This algorithm will find the optimal rotation to minimize the absolute mean difference between the maps.
# It works by computing the absolute mean difference between two evaluated interpolation grids (across a set of rotation values, user defined by *theta*).
# This technique must be carefully monitored to make sure the orientation chooses a reasonable fix.
# To visualize rotation optimization, we've provided an easy interface for the rotation_optimize function, which can be invoked with the sample code:
RESOLUTION_var <- 100
spin_spatial_res <- 32
theta_var <- 1
rotation_df <- rotation_optimize(Pmol_753,
								 Pmol_752,
								 spatial_res=spin_spatial_res,
								 theta_interval=theta_var
								)
plot_rotation_optimize(rotation_df)
rotate_op <- optimal_rotation(rotation_df)[1]
