RESOLUTION_var = 265
my_lambda = 0.01
retinaplot(Pmol_752, spatial_res=RESOLUTION_var, polynomial_m = 3) #quick big
retinaplot( Pmol_752,
			spatial_res=RESOLUTION_var,
			contour_breaks_source  =  c(0,100000), 
			col_breaks_source      =  c(0,100000),
			col_levels=50,
			contour_levels=20,
			rotation=0,
			inner_eye_view=TRUE,
			lambda=.001,
			polynomial_m=NULL)

retinaplot( Pmol_752,
			spatial_res=RESOLUTION_var,
			contour_breaks_source  =  c(0,100000), 
			col_breaks_source      =  c(0,100000),
			col_levels=50,
			contour_levels=20,
			rotation=0,
			inner_eye_view=TRUE)




#Change lambda and the polynomial


retinaplot( Pmol_752,
			spatial_res=RESOLUTION_var,
			contour_breaks_source  =  c(0,100000), 
			col_breaks_source      =  c(0,100000),
			rotation=-30,
			inner_eye_view=TRUE)

# retinaplot(Pmol_753,contour_breaks_source  =  c(0,100000), 
# 					col_breaks_source      =  c(0,100000), inner_eye_view=TRUE)
# retinaplot(Ntae_381,contour_breaks_source  =  c(0,100000), 
# 					col_breaks_source      =  c(0,100000), inner_eye_view=TRUE)