

context('test_smoothing_params')

test_that("we can compute different versions of the retinaplot smoothing with lamdba and polynomial values", {

	

	RES_variable <- 512

pdf(file="smoothing_combos.pdf", width=8, height=5.5, family="Palatino")
	poly <-c(2, 3, 4, 5)
	lambda <- c(1e-2,1e-3,1e-4,1e-5)
	# Generate error vectors for each combination
	# Save each plot to a PDF in the working directory

	source("3hgbqg" %>% path_to_main_file_for_test_retina)
retina_A <- main_3hgbqg("3hgbqg" %>% get_path_to_test_retina_folder)


ERR_dat <- polynomial_vs_lambda(retina_A, spatial_res= RES_variable, polynomial_m_vec=poly, lambda_vec=lambda)

sampling_location_number = length(retina_A$fit_data1$y)

ERR_dat <- ERR_dat[-1] #take off the NA placeholder column

X_median	<- retina::reorder_columns(ERR_dat, median)
X_min 		<- retina::reorder_columns(ERR_dat, min)
X_max  		<- retina::reorder_columns(ERR_dat, max)
X_mean  	<- retina::reorder_columns(ERR_dat, mean)
X_sd  		<- retina::reorder_columns(ERR_dat, sd)
X_range_len <- retina::reorder_columns(ERR_dat, retina::range_len)

the_types <- c( 'X_median',
				'X_min',
				'X_max',
				'X_mean',
				'X_sd',
				'X_range_len')

for (e in the_types) {
	err_obj <- eval(parse(text=e))
	boxplot(err_obj, horizontal=TRUE, las=1, col='lightgray', cex=.2, pch=20, sub=paste("n = ", sampling_location_number), cex.axis=0.4, main=paste0("Ordered by ", e))
}

boxplot(X_sd, horizontal=TRUE, las=1, col='lightgray', asp=1, cex=.2, pch=20, sub=paste("n = ", sampling_location_number), cex.axis=0.4, main=paste0("Ordered by sd"))

barplot(apply(X_sd, 2 ,sd), horiz=TRUE, las=1)
dev.off()

})