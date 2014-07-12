##alpha.

retina_list = list(Pmol_753=Pmol_753, Pmol_752=Pmol_752, Ntae_381=Ntae_381)
ret_list_len <- length(retina_list)
total_entries <- ret_list_len^2
mat_vec <- c()
counter=1
for (map1 in names(retina_list)) {
	for (map2 in names(retina_list)) {
		
		message(paste(counter, 'of',total_entries))
		rotation_df <- rotation_optimize(retina_list[[map1]], 
										 retina_list[[map2]], 
										 spatial_res    = 4,
										 theta_interval = 10)
		map1_map2_rotation_value<- as.numeric(optimal_rotation(rotation_df)[1])
		message(paste(map2 , 'is',
					  map1_map2_rotation_value,'degrees off, when',
					  map1, "is fixed" ))
	    print(map1_map2_rotation_value)
	    mat_vec <- c(mat_vec, map1_map2_rotation_value)
	    counter=counter+1
	}

}

MAT <- t(matrix(mat_vec, ret_list_len, ret_list_len)) #coerce the vector into matrix


partial_av<- function(retina_object, rotation, spatial_res, reflect=TRUE, ...){
	x2 <- retina_object$azimuthal_data.datapoints[[1]]$x
	y2 <- retina_object$azimuthal_data.datapoints[[1]]$y
	z2 <- retina_object$azimuthal_data.datapoints[[1]]$z

	# Rotation
	ori <- c(-90, 90, rotation)
	az <- mapproject(x=retina_object$trimmed_data[,2], retina_object$trimmed_data[,1], 
					   projection="azequidistant", orientation=ori)
	x2 <- az$x #overwirte old value
	ifelse(reflect,x2<-(-1.0*x2), message(''))
	y2 <- az$y #overwrite old value

	# Fitting models
	map2fit <- PolarImageInterpolate(x2,y2,z2, spatial_res,
									plot_suppress=TRUE, extrapolate=TRUE,
									outer.radius=pi/2.0,
									falciform_coords = retina_object$azimuthal_data.falciform[[1]],...)
	MAT <- map2fit[[2]]$z	
	#reflect the eye
	return(MAT/n)
}
partial_Pmol_752 <- partial_av(Pmol_752, 1, spatial_res=100)
partial_Pmol_753 <- partial_av(Pmol_753, 1, spatial_res=100)
partial_Ntae_381 <- partial_av(Ntae_381, 1, spatial_res=100)
retinaplot(Pmol_752)


map_vec_average <- function(retina_list, rotation_indices){
	sum_mat <- partial_av(retina_list[[1]], n, spatial_res=1000)
	n <- length(retina_list)
	for (i in 2:n) {
		sum_mat <- sum_mat + partial_av(retina_list[[i]], n, spatial_res=1000)
		browser()
	}
	return(sum_mat)
}

sum <- map_vec_average(list(Pmol_752,Pmol_753, Ntae_381), c(0,0,179))

plot_from_MAT(MATRIX=sum)
