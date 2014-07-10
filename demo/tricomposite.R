##alpha.

retina_list = list(Pmol_753, Pmol_752, Ntae_381)

	for (map1 in names(retina_list)) {
	for (map2 in names(retina_list)) {
		rotation_df <- rotation_optimize(retina_list[[map1]], 
										 retina_list[[map2]], 
										 spatial_res    = 32,
										 theta_interval = 30)
		map1_map2_rotation_value <- optimal_rotation(rotation_df)
		paste(map2 , 'is',map1_map2_rotation_value,'degrees off, when' ,map1, "is fixed" )
	    print(map1_map2_rotation_value)
	}
}
