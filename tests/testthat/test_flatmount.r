context('flatmount for comparison')
test_that('flatmounts can be produced for comparison and inclusion in supplementary material', {
    path_to_retina_data_folder <- "/Users/briancohn/Documents/GitHub/bc/retina/inst/extdata/asalmon_mcw/diagram_retina"
	retistruct_dataset <- retistruct.read.dataset(path_to_retina_data_folder)
    retistruct.reconstruct(retistruct_dataset)
    flatplot(retistruct_dataset, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=TRUE, datapoints=TRUE, grid=FALSE, landmarks=TRUE)
    
}	