# Retina Package Tutorial
=====


## Install:
[README Instructions](README.md "Readme instructions on bc/retina")  
You should be able to run this successfully:
```R
test()
```

## Questions, Errors, and Comments?
[Submit a quick help ticket](https://github.com/bcohn12/retina/issues/new "") or brian.cohn@usc.edu


# Retinal data pre-processing
=====

##1. Manually set up a folder containing your starting retina files
###Create a folder called 'diagram_retina', which will contain at first:  
1. `diagram_retina_screenshot.png` - A screen shot of the stereology software, showing the sampling locations on top of the retinal outline. [Example PNG](https://github.com/bcohn12/retina/blob/master/inst/extdata/tutorial_data/diagram_retina_screenshot.png "")

2. `site_counts_from_stereology.csv` - a 2 column file with cell counts for each sampling location. This comes from your stereology software or from manual recording. [Example Site Count CSV](https://github.com/bcohn12/retina/blob/master/inst/extdata/tutorial_data/site_counts_from_stereology.csv "")

##2. Manually make `xyz.csv`, a file of cell counts with their unique `x` and `y` positions.
You’ll do this by using the `site_counts_from_stereology.csv` and the `diagram_retina_screenshot.png` to assign `x` and `y` values.

![Image of mapping to xyz](tutorial_pix/make_xyz_mapping.png "Example")
>**Example 1 (above)** Sampling location 13 is at the location `x=5`, `y=4`, because the sampling location is in the top left center of the box in the 4th row, at the 1st column. The `site_counts_from_stereology.csv` shows us that `z=10`.  

>**Example 2** Sampling location 5 is at (7,5). The sampling file says there's 3 cells at location 5. The row in xyz.csv for this sampling location would be `x=7`,`y=5`,`z=3`.  

The origin `(1,1)` is the bottom left corner of the grid, where `x` is the column (vertical), and `y` is the row (horizontal).  
You should finish filling out a three column csv, with the three columns named (without quotes) "x" "y" "z" . Number of rows (excluding the header) will be equal to the number of sampling locations.  
![Image of xyz.csv_example](tutorial_pix/xyz_template.png "Example")

Save this as `xyz.csv` in the 'diagram_retina' folder.

##3. Use ImageJ to record outline vertices and calibrate datapoints.
Install [ImageJ](http://imagej.nih.gov/ij/download.html "Download link") (for PC/Mac/Linux.  
Open `diagram_retina_screenshot.png`

###Trace the retinal outline to make an outline.ROI file
The file name has to be outline.roi.
![Select the polygon](tutorial_pix/select_polygon.png "Select the polygon tool in ImageJ")  
Select the polygon selection tool  
[Select](tutorial_pix/select_polygon.png) the Polyton Selection tool and [trace the polygon.](tutorial_pix/polygon_tracing.png) Finish the polygon by clicking on the first point. Be precise with the incision slits.

From the menubar `Analyze>Tools>ROI Manager`. Press the `Add [t]` [button](https://cloud.githubusercontent.com/assets/4623063/12045123/027fcb7a-ae52-11e5-99c0-42de86a595e3.gif).
Press the `More` button, then `Save` it as `outline.roi` in your `/diagram_retina` folder.

[Trace the falciform process/optic disk with the polygon tool](tutorial_pix/falc_traced.png) and menubar `File > Save As > Save as XY Coordinates`; save as `falc.txt` in diagram_retina.

# If your falciform coordinate values are totally different than the values of the sampling site coordinates:
Via the menubar: `Analyze > Set Scale` and click to **Remove Scale**. This will make sure that the values in `falc.txt` are in pixel units.

### Record the ImageJ coordinates of the sampling location bounds
![Minx](tutorial_pix/minx.png "")  
Find the coordinates of the outermost sampling locations. Hover with your mouse and look at the live-updated (x,y) coordinates in the tool bar.  
**minX Example** Sampling location 10 is the furthest point to the left, and it's at `x=25`.  

**maxX Example** Sampling locations 5, 15, and 16 are on the far right, and share the same edge. They all have `x=566`, so set `maxX` to `566`.  

**deltaX Example** deltaX=(maxX-minX)/6, because there are six vertical lines.

**deltaY Example** deltaY=(maxY-minY)/5, because there are five horizontal lines.

Delta X and Delta Y are positive numbers.

**NOTE** Convert both minY and maxY so they are negative numbers (even though they may look like they're positive in ImageJ.
Write down the minX, maxX, minY, maxY somewhere safe. You'll use it when you call the retina_object function.




4. Markup the tears
=====
View the indices of the outline coordinates. If you made an outline of 25 points, then the indices will be between 1 and 25. This is how we tell R where to stitch up the wholemount incision tears. It's like defining where to stitch.
```r
outline_coordinates <- tear_markup_plot(path_to_retina_data_folder)
```

```r
#V0 point of the tear, the point before it (smaller than the V point), point after it]
#"V0","VB","VF"
#these are the actual points for Pmol_752, as an example. You can practice on that retina in `inst/extData/Pmol_752`
tear_coordinates_dataframe <- rbind(
	c(31,33,30),
	c(2,3,1),
	c(22,23,21),
	c(15,17,14),
	c(10,11,9)
																	)
assemble_tear_file(tear_coordinates_dataframe, path_to_retina_data_folder)
assemble_point_coordinates_file(outline_coordinates, path_to_retina_data_folder)
assemble_markup_file('left', path_to_retina_data_folder, nasal_outline_index=27, dorsal_outline_index=NA, phi0=retinal_phi0(ED=1.0, AL=0.95))
```

### Optional: Adding retinal perimeter latitude for highly non-hemispherical retinae
If you would like to set the Retinal perimeter latitude, paste it into the `phi0` box (in units of degrees).
We include the `retinal_phi0` function to help users make an approximation of this value with recorded eye measurements. Contact brian.cohn@usc.edu if you want this feature enabled.

#### Sample calculation
```r
> ED = 4.8
> AL = 3.43
> phi0 <- retinal_phi0(ED, AL)
> phi0
> [1] 25.41469
```

Set the working directory to the folder which contains the *diagram_retina* folder.

Set your parameters
```R
LAMBDA_var <- 0.001 #Thin plate spline smoothing. Lambda=0 would be interpolation.
RESOLUTION_var <- 500 #desired plot width in pixels
IJ <-data.frame(maxX = 1437,  # rightmost counting location's X value
				maxY = -45,   # topmost counting location's Y value
				minX = 469,   # leftmost counting location's X value
				minY = -1029, # bottommost counting location's Y value
				deltaX = 60,  # average ImageJ pixel distance (in the X axis) between counting locations.
				deltaY = 61)  # average ImageJ pixel distance (in the Y axis) between counting locations.
```

Eye measurements from dissection are not required for retinal mapping, so if these were not collected or are unavailable, simply put 1.0 for ED, AL, LD.
```R
my_retina <- retina_object(
	path = "diagram_retina",

	#Eye Measurements from dissection
		ED = 4.8,    #Eye diameter (mm)
		AL = 3.43,   #Eye axial length (mm)
		LD = 1.8,    #Eye lens diameter (mm)

	#Stereology Parameters
		height = 25 ,# height of the counting frame in microns
		width  = 25, # width of the counting frame in microns

	#Plotting Parameters
		lambda = LAMBDA_var,          #see fields::Tps for more information
		extrapolate = TRUE ,          #Predicts densities to the equator.
		spatial_res = RESOLUTION_var,
		rotation_ccw = -90, # when set to -90 degrees, the rotation is unaltered from the measured orientation.
		plot_suppress=TRUE,
	#ImageJ Datapoint Calibration Measurements
		IJcoords = IJ)
```
# Visualization and Diagnostics
```R
retinaplot(my_retina) ##Plot the retina
retinaplot(my_retina, spatial_res = 1000) ##Plot a 1000x1000 pixel sized-retinaplot (long computation time, high resolution)
retinaplot(my_retina, spatial_res = 100) ##Plot a 100x100 pixel sized-retinaplot (short computation time, low resolution)
#try both of these to make sure you aren't predicting too much to the equator.
retinaplot(my_retina, extrapolate=FALSE)
retinaplot(my_retina, extrapolate=TRUE)
the retina that was sampled.
retinaplot(my_retina, contour_levels=20) ##Define how many topographic lines you want to define the contours.
These are fitting parameters you can also modify:
`lambda`
`polynomial_m`

# Compute the retinal perimeter, to define the width of the plot in mm.
ret_perimeter_len <- semi_ellipse_perimeter(a=ED, b=AL)/2

# Fit diagnostics
fit_plots(my_retina$fit_data1)
```
w
# Further functionality
```R
?composite_map 				#make an average of 2 maps
?polynomial_vs_lambda 		#compare smoothing parameters on your map
```

Users can plot multiple falciform processes atop one another by extracting them from the retinal composite, and plotting the polygons:
```R
# where fc1 is the first retina’s falciform/optic disc coordinates
polygon(fc1[,1], fc1[,2], col=rgb(0, 0, 0,0.5), lty="solid", border="gray42")
polygon(fc2[,1], fc2[,2], col=rgb(0, 0, 0,0.5), lty="solid", border="gray42")
```

# Saving to a paginated PDF
```R
pdf("maps.pdf", width=8.5, height=6)
	#Plot your retinal map(s)
dev.off()
```

# Bugs people have encountered and how to fix them:
1. `Error in Summary.factor(c(34L, 15L, 240L, 225L, 239L, 14L, 33L, 101L,  : ‘min’ not meaningful for factors` . Make sure your `xyz.csv` file is 'comma separated', and not colon-separated.
