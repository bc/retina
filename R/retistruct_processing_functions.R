##' @title Retistruct Wrapper
##' @description
##' Reads in the folder, reconstructs, and returns the retinal object.
##' @param path string path to the folder with retinal data.
##' @return rad Retinal data in radian units.
##' @author Brian Cohn \email{brian.cohn@@usc.edu}
##' @references Sterratt et. al. 2013
dss_retistruct_processing <- function(path) {
    do_not_print <- function(string){}
    rad <- retistruct.read.markup(retistruct.read.dataset(path))
    rad <- retistruct.reconstruct(rad, report=do_not_print, plot.3d=FALSE, )  ## Reconstruct (computation intensive)
    return(rad)
}

##' @title Set up tear matrix
##' @description
##' takes in a list of tear triplets (numeric vectors of 3 elements) and outputs the dataframe, ready for retistruct processing
##' @return mat matrix of tear vals, ready for T.csv saving
##' @author Brian Cohn \email{brian.cohn@@usc.edu}
##' @export
compose_tear_triplets_dataframe <- function(list_of_tear_triplets) {
	permute_tear_vertices <- function(list_of_tear_triplets){
		lapply(list_of_tear_triplets, function(x){
		retistruct::labelTearPoints(annotated_outline_with_tears,x)
			})
		}

    	do.call('rbind',permute_tear_vertices(list_of_tear_triplets))
}
