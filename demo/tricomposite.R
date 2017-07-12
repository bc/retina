
retina_list <- list(Pmol_752=Pmol_752,
					Pmol_753=Pmol_753,
					Pmol_751=Pmol_752 )

ret_list <- vector_retina_composite(retina_list,
									plot=TRUE,
									spin_spatial_res=8,
								  	plot_spatial_res=8,
								  	theta_interval=180,
                    rotation_degree_list=as.matrix(c(0,0,0)))

ret_list <- vector_retina_composite(retina_list,
									plot=TRUE,
									spin_spatial_res=8,
								  	plot_spatial_res=8,
								  	theta_interval=180,
                    rotation_degree_list=FALSE)
