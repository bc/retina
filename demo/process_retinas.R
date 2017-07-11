#Retinal Plotter Example Code: Brian Cohn March 28, 2014
#Bring in package/function files
require(testthat);
require(sphereplot); require(mapproj); require(rgl);
require(fields); require(RColorBrewer);
require(retina)

#variable numbers
RESOLUTION_var <- 512
SPIN_resolution <- 64
my_lambda <- 0.001

# Run retistruct()
# Markup the retina's incisions
# Set eye to Left or Right eye
# Set a most dorsal, or most nasal point on the retina
# 'Save' and close the retistruct GUI
Ntae_381_coords <- data.frame(  minX   = 42,  #the leftmost counting frame's X value
								maxX   = 597, #the rightmost counting frame's X value (from the ImageJ outline PNG)
								deltaX = ((597-42)/17), #the average ImageJ pixel distance (in the X axis) between counting frame locations from the outline image.
								minY   = -584, #bottommost counting frame's Y value
								maxY   = -32, #the topmost counting frame's Y value
								deltaY = (584-32)/17   ) #the average ImageJ pixel distance (in the Y axis) between counting frame locations from the outline image.

Ntae_381<- retina_object(
	path = system.file(package = "retina", 'extdata/Ntae_381'),
	#Eye Measurements
	 	LD = NULL,
		ED = 4.508 , #ED EYE DIAMETER (at the hemisphere) in mm
		AL = 3.1395 ,#AL Axial Length (from front to back of eye) in mm
	#Stereology Parameters
		height = 25 ,# height of the counting frame in microns
		width = 25, # width of the counting frame in microns
	#Plotting Parameters
		lambda=my_lambda ,#defines the smoothing parameter for thin plate spline interpolation (see fields::Tps for more information)
		extrapolate=TRUE ,#when true, the plotter will create a full representation of the hemisphere even outside of the bounds of measured points.
		spatial_res=RESOLUTION_var, #number of pixels wide the plotter will use. a value of 1000 will create a 1000x1000 pixel plot
		rotation_ccw= -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
	#ImageJ Datapoint Calibration Measurements
	IJcoords = Ntae_381_coords,
	plot_suppress=TRUE
	#End parameter Specification
)
Pmol_753_coords <- data.frame(  maxX=1437,
								maxY=-45,
								minX=469,
								minY=-1029,
								deltaX=60,
								deltaY=61) #the average ImageJ pixel distance (in the Y axis) between counting frame locations from the outline image.
Pmol_753<- retina_object(
	path = system.file(package = "retina", 'extdata/Pmol_753'),
	#Eye Measurements
	ED = 4.8,
	AL = 3.43,
	LD = 1.8,
	#Stereology Parameters
		height = 25 ,# height of the counting frame in microns
		width  = 25, # width of the counting frame in microns
	#Plotting Parameters
		lambda = my_lambda ,#defines the smoothing parameter for thin plate spline interpolation (see fields::Tps for more information)
		extrapolate = TRUE ,#when true, the plotter will create a full representation of the hemisphere even outside of the bounds of measured points.
		spatial_res = RESOLUTION_var, #number of pixels wide the plotter will use. a value of 1000 will create a 1000x1000 pixel plot
		rotation_ccw = -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
		plot_suppress=TRUE,
	#ImageJ Datapoint Calibration Measurements
	IJcoords = Pmol_753_coords)




Pmol_752_coords <- data.frame(maxX=848,
								maxY=-48,
								minX=40,
								minY=-965,
								deltaX=(848-40)/15,
								deltaY=(965-48)/17) #the average ImageJ pixel distance (in the Y axis) between counting frame locations from the outline image.

## this part demonstrates how to combine sampling site count data, with your manually recorded XY location of each sampling site.
	ssites<- read.csv(system.file(package = "retina", 'extdata/Pmol_752/Pmol_752_ssite_locations.csv'))
	counts<- read.csv(system.file(package = "retina", 'extdata/Pmol_752/Pmol_site_counts.csv'))
	colnames(ssites) <- c("x", "y", "samplingsite")
	colnames(counts) <- c("samplingsite", "count")
	Pmol_752_xyz <- ssite_merge(ssites,counts)
	colnames(Pmol_752_xyz) <- c("x","y","z")
	write.csv(Pmol_752_xyz, system.file(package = "retina", 'extdata/Pmol_752/xyz.csv'), row.names=FALSE)
Pmol_752<- retina_object(
	path = system.file(package = "retina", 'extdata/Pmol_752'),
	ED=5.225,
	AL=3.9,
	LD=2.05,
	height = 25,
	width = 25,
	lambda = my_lambda,
	extrapolate=TRUE,
	spatial_res = RESOLUTION_var,
	rotation_ccw = -90,
	plot_suppress=TRUE,
	IJcoords = Pmol_752_coords)
