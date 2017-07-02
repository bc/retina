##' @title Compute rotation matrix
##' @description
##' Fixes each retina in the inputted list of retina-objects, generating the optimal rotation (in degrees) with respect to the other retinas.
##' @param retina_list list of retinal objects. for example, list(Pmol_753=Pmol_753, Pmol_752=Pmol_752, Ntae_381=Ntae_381)
##' @param spatial_res width of the evaluation grid that's used for every combination. 
##' @param theta_interval This value (an integer in degrees) will generate (360/theta_interval) rotations to evaluate.
##' @return rotation_mat Matrix where columns are fixed, rows are rotated. Column 1 is the first element in the retina_list
##' @author Brian Cohn \email{brian.cohn@@usc.edu} and Lars Schmitz
##' @family spin_optimization
##' @export
compute_rotation_matrix <- function(retina_list, spatial_res = 16, theta_interval = 10) {
    ret_list_len <- length(retina_list)
    total_entries <- ret_list_len^2
    mat_vec <- c()
    counter = 1
    for (map1 in names(retina_list)) {
        for (map2 in names(retina_list)) {
            message(paste(counter, "of", total_entries))
            rotation_df <- rotation_optimize(retina_list[[map1]], retina_list[[map2]], 
                spatial_res = 4, theta_interval = 10)
            map1_map2_rotation_value <- as.numeric(optimal_rotation(rotation_df)[1])
            message(paste("when", map1, "is fixed", map2, "is", map1_map2_rotation_value, 
                "degrees off."))
            mat_vec <- c(mat_vec, map1_map2_rotation_value)
            counter = counter + 1
        }
    }
    
    MAT <- t(matrix(mat_vec, ret_list_len, ret_list_len))  #coerce the vector into matrix
    return(MAT)
}
