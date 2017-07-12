
RESOLUTION_var = 512
my_lambda = 0.001
retinaplot(Pmol_752, spatial_res=RESOLUTION_var, polynomial_m = 2) #quick big
retinaplot( Pmol_752,
			spatial_res=RESOLUTION_var,
			rotation=0,
			inner_eye_view=TRUE,
			lambda=.001,
			polynomial_m=2)

retinaplot( Pmol_752,
			spatial_res=RESOLUTION_var,
			rotation=0,
			inner_eye_view=TRUE)




#Change lambda and the polynomial


retinaplot( Pmol_752,
			spatial_res=RESOLUTION_var,
			rotation=30,
			inner_eye_view=TRUE)
