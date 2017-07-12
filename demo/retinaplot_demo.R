par(mfrow=c(1,1))
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

# retinaplot(Pmol_753,contour_breaks_source  =  c(0,100000),
# 					col_breaks_source      =  c(0,100000), inner_eye_view=TRUE)
# retinaplot(Ntae_381,contour_breaks_source  =  c(0,100000),
# 					col_breaks_source      =  c(0,100000), inner_eye_view=TRUE)
