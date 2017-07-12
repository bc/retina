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
