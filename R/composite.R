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
