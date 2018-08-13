

context("test_spin_optimization")

test_that("test_spin optimization", {
    # Lets make a species average.  Use the *composite_map* function to make a
    # species average of two retinae.  The user can enable 'spin optimization', which
    # we designed to help users overlay maps with minimal orientation error.  This
    # algorithm will find the optimal rotation to minimize the absolute mean
    # difference between the maps.  It works by computing the absolute mean
    # difference between two evaluated interpolation grids (across a set of rotation
    # values, user defined by *theta*).  This technique must be carefully monitored
    # to make sure the orientation chooses a reasonable fix.  To visualize rotation
    # optimization, we've provided an easy interface for the rotation_optimize
    # function, which can be invoked with the sample code:
    RESOLUTION_var <- 100
    spin_spatial_res <- 32
    theta_var <- 1
    
    retina_A <- main_3hgbqg("3hgbqg" %>% get_path_to_test_retina_folder)
    retina_B <- main_40oik5("40oik5" %>% get_path_to_test_retina_folder)
    
    rotation_df <- rotation_optimize(retina_A, retina_B, spatial_res = spin_spatial_res, 
        theta_interval = theta_var)
    plot_rotation_optimize(rotation_df)
    rotate_op <- optimal_rotation(rotation_df)[1]    
})