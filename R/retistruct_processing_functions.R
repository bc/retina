##' @title Retistruct Wrapper
##' @description
##' Reads in the folder, reconstructs, and returns the retinal object.
##' @param path string path to the folder with retinal data.
##' @return rad Retinal data in radian units.
##' @author Brian Cohn \email{brian.cohn@@usc.edu}
##' @references Sterratt et. al. 2013
##' @importFrom retistruct retistruct.read.markup retistruct.reconstruct retistruct.read.dataset
##' @export
dss_retistruct_processing <- function(path) {
    do_not_print <- function(string) {
    }
    rad <- retistruct.read.markup(retistruct.read.dataset(path))
    rad <- retistruct.reconstruct(rad, report = do_not_print, plot.3d = FALSE)  ## Reconstruct (computation intensive)
    return(rad)
}

##' @title Set up tear matrix
##' @description
##' takes in a list of tear triplets (numeric vectors of 3 elements) and outputs the dataframe, ready for retistruct processing
##' @param list_of_tear_triplets list of 3-element vectors
##' @param outline_object See ?retistruct::labelTearPoints for documentation
##' @return mat matrix of tear vals, ready for T.csv saving
##' @author Brian Cohn \email{brian.cohn@@usc.edu}
##' @export
##' @importFrom retistruct labelTearPoints
compose_tear_triplets_dataframe <- function(list_of_tear_triplets, outline_object) {
    permute_tear_vertices <- function(list_of_tear_triplets, outline_object) {
        my_outline <- outline_object
        lapply(list_of_tear_triplets, function(x) {
            labelTearPoints(my_outline, x)
        })
    }
    do.call("rbind", permute_tear_vertices(list_of_tear_triplets, outline_object))
}
##' @title Set up tear matrix within the annotated outline object
##' @description
##' assembles the outline object with tear data at $V0, $VB, and $VF
##' @param outline_object See ?retistruct::labelTearPoints for documentation
##' @param tear_coordinates_dataframe dataframe with $V0, $VB, and $VF columns as produced by compose_tear_triplets_dataframe
##' @return annotated_outline See ?retistruct::labelTearPoints for documentation
##' @author Brian Cohn \email{brian.cohn@@usc.edu}
##' @export
update_outline_object_tears <- function(outline_object, tear_coordinates_dataframe) {
    outline_copy <- outline_object
    outline_copy$V0 <- tear_coordinates_dataframe[, 1]
    outline_copy$VB <- tear_coordinates_dataframe[, 2]
    outline_copy$VF <- tear_coordinates_dataframe[, 3]
    return(outline_copy)
}

##' @title Generate outline with tears
##' @description
##' Creates a valid AnnotatedOutline, ready for retistruct
##' @param path_to_retina_data_folder path to working directory
##' @param list_of_tear_triplets list of 3-element vectors
##' @param outline_coordinates Outline coordinates XY str generated from tear_markup_plot function
##' @return outline_with_tears AnnotatedOutline. See ?retistruct::AnnotatedOutline
##' @importFrom retistruct AnnotatedOutline Outline
##' @author Brian Cohn \email{brian.cohn@@usc.edu}
##' @export
generate_outline_with_tears <- function(outline_coordinates, list_of_tear_triplets, 
    path_to_retina_data_folder) {
    unannotated_outline <- Outline(outline_coordinates, scale = NA, im = NULL)
    outline_object <- AnnotatedOutline(unannotated_outline)
    outline_with_tears <- update_outline_object_tears(outline_object, assemble_tear_file(compose_tear_triplets_dataframe(list_of_tear_triplets, 
        outline_object), path_to_retina_data_folder))
    return(outline_with_tears)
}
