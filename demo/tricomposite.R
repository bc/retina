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
		map1_map2_rotation_value <- as.numeric(optimal_rotation(rotation_df)[1])
		message(paste(map2 , 'is',map1_map2_rotation_value,'degrees off, when' ,map1, "is fixed" ))
	    print(map1_map2_rotation_value)
	    mat_vec <- c(mat_vec, map1_map2_rotation_value)
	    counter=counter+1
	}

}


MAT <- t(matrix(mat_vec, ret_list_len, ret_list_len))


vector_maps_composite <- function(retina_list, MAT){
	rotations <- MAT[1,]

}


rotate_retina_object <- function(retina_object, rotation){
	
}