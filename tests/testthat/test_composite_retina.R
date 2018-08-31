
context("test_composite_retina")

test_that("we can combine 5 retinas into one", {
    spatial_res <- 10
    retina_A <- main_3hgbqg("3hgbqg" %>% get_path_to_test_retina_folder %>% map_to_matrix(spatial_res, projection_type='azequidistant'))
    retina_B <- main_40oik5("40oik5" %>% get_path_to_test_retina_folder %>% map_to_matrix(spatial_res, projection_type='azequidistant'))
    retina_C <- main_8pm223("8pm223" %>% get_path_to_test_retina_folder %>% map_to_matrix(spatial_res, projection_type='azequidistant'))
    retina_D <- main_gv3igs("gv3igs" %>% get_path_to_test_retina_folder %>% map_to_matrix(spatial_res, projection_type='azequidistant'))
    retina_E <- main_jhvbwt("jhvbwt" %>% get_path_to_test_retina_folder %>% map_to_matrix(spatial_res, projection_type='azequidistant'))
    
    # Example Code:
    list_of_retina_matrices <- list(retina_A,
                                    retina_B,
                                    retina_C,
                                    retina_D,
                                    retina_E)
    
    mean_map_mat <- multimap_composite(list_of_retina_matrices, map_names_vector=c("3hgbqg","40oik5","8pm223","gv3igs","jhvbwt"))
    plot_multimap_composite(mean_map_mat)

})

test_that("we can combine two retinas into one", {
    # If spin_optimization is disabled, no rotation will be performed.  composite_map
    # evaluates the density across an interpolation grid at the desired *spatial_res
    # = n* spatial resolution, and saves this data into an n by n matrtix.  Finally,
    # both matrices are added, and divided by two to produce the 'matrix average'.
    # This map can be plotted by passing the matrix to the plot_from_MAT function.
    # **note: the resultant map does not contain a combined falciform process, and
    # the original sampling site points are not visible on the species average.
    RESOLUTION_var <- 512
    retina_A <- main_3hgbqg("3hgbqg" %>% get_path_to_test_retina_folder)
    retina_B <- main_40oik5("40oik5" %>% get_path_to_test_retina_folder)
    
    
    retina_composite <- composite_map(retina_A, retina_B, rotation = TRUE, spatial_res = RESOLUTION_var)
    
    # retina_composite <- composite_map(retina_B, retina_A, rotation=TRUE,
    # spatial_res=RESOLUTION_var)
    boxplot(retina_A$azimuthal_data.datapoints[[1]]$z, retina_composite, retina_B$azimuthal_data.datapoints[[1]]$z, 
        names = c("Sample 752", "Average", "Sample 753"), xlab = "Retinal ganglion cells per square mm",
        horizontal = TRUE, pch = 20, cex = 0.5, col = "lightgray")
    rho_max <- max(retina_composite)
    rho_min <- min(retina_composite)
    # plot the species composite
    retinaplot(retina_A, contour_breaks_source = c(rho_min, rho_max), col_breaks_source = c(rho_min, 
        rho_max), col_levels = 50, contour_levels = 20, spatial_res = RESOLUTION_var)
    retinaplot(retina_B, contour_breaks_source = c(rho_min, rho_max), col_breaks_source = c(rho_min, 
        rho_max), col_levels = 50, contour_levels = 20, spatial_res = RESOLUTION_var)
    retinaplot(retina_B, contour_breaks_source = c(rho_min, rho_max), col_breaks_source = c(rho_min, 
        rho_max), col_levels = 50, contour_levels = 20, spatial_res = RESOLUTION_var, 
        rotation = -30)
    plot_from_MAT(retina_composite, extrapolate = FALSE, spatial_res = RESOLUTION_var, 
        col_levels = 50, contour_levels = 20, contour_breaks_source = c(rho_min, 
            rho_max), col_breaks_source = c(rho_min, rho_max), z1 = retina_B$azimuthal_data.datapoints[[1]]$z, 
        z2 = retina_A$azimuthal_data.datapoints[[1]]$z)
})



