##' @title Generate a retinal composite
##' @description Combine multiple retinal maps together. This computes a thin plate spline for each of the retinas, decides on a common XY grid to compare points at, and predicts each of the maps at those locations. Then, it takes the average of the density values at those locations.
##' @param retina_list A list containing retina objects, with names(retina_list) all defined.
##' @param spin_spatial_res Resolution passed to the spin-optimization.
##' @param theta_interval number of degrees to traverse in spin-optimization
##' @param plot_spatial_res Resolution passed to the final plotter
##' @param plot Logical, by default false, but when true it will output the plot where the first element in the list is fixed.
##' @param ... Arguments passed from other functions. You can access the thin plate spline interpolator here.
##' @return composite_mat Matrix of the composite retina, where the first element of the retina list is fixed.
##' @author Brian Cohn \email{brian.cohn@@usc.edu} and Lars Schmitz
##' @export
vector_retina_composite <- function(retina_list, spin_spatial_res = 128, plot_spatial_res = 512,
    theta_interval = 10, plot = FALSE, rotation_degree_list, ...) {
    if (class(rotation_degree_list) != "matrix") {
        rotation_degree_list <- compute_rotation_matrix(retina_list, spatial_res = spin_spatial_res,
            theta_interval = theta_interval, ...)[1, ]
        message("No coordinates submitted, spin optimization chosen")
    }
    n <- length(retina_list)
    sum <- map_vec_sum(retina_list, rotation_degree_list, spatial_res = plot_spatial_res)
    composite_mat <- sum/n  #get the average
    if (plot)
        {
            message("Plotting")
            dev.new()
            plot_from_MAT(MATRIX = composite_mat, extrapolate = TRUE, spatial_res = plot_spatial_res,
                col_levels = 50, contour_levels = 20, contour_breaks_source = c(min(composite_mat),
                  max(composite_mat)), col_breaks_source = c(min(composite_mat),
                  max(composite_mat)), z1 = Pmol_753$azimuthal_data.datapoints[[1]]$z,
                z2 = Pmol_752$azimuthal_data.datapoints[[1]]$z)
        }  #end plotif
    return(composite_mat)
}

##' @title Compute the matrix sum of multiple retinas
##' @description This is a subfunction, that is useful when creating average maps.
##' @param retina_list A list containing retina objects, with names(retina_list) all defined.
##' @param rotation_indices Each value in this vector will correspond to the amount of rotation applied to each retina in the list.
##' @param ... Arguments passed from other functions. You can access the thin plate spline interpolator here.
##' @return sum_mat The sum of all retina objects in the list, with respect to the rotation indices.
##' @author Brian Cohn \email{brian.cohn@@usc.edu} and Lars Schmitz
map_vec_sum <- function(retina_list, rotation_indices, ...) {
    sum_mat <- mat_from_ret_obj(retina_list[[1]], rotation = 0, ...)
    n <- length(retina_list)
    for (i in 2:n) {
        sum_mat <- sum_mat + mat_from_ret_obj(retina_list[[i]], rotation_indices[i],
            ...)
    }
    return(sum_mat)
}
