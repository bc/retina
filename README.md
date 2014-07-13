retina in R
======

Open-source tools for visualizing and comparing retinal cell density data.


> A _retinal ganglion cell map_ helps biologists visualize the distribution of receptive cells across the rear surface of the vertebrate eye. This topographic map shows density in a way which highlights areas of higher visual performance.

Make a retinal map:
======
After surveying cell density across the retina using a stereology-equipped microscope:

Use _retina_ to

      + Import retinal outline and sampling data
      + Reinstate the hemispherical shape of the data (1)
      + Make species average maps
      + Evaluate parameters and fit of smoothing models
      + Plot a contoured topographic heatmap

Plots were carefully designed to be readable by viewers who may have colorblindness.


Install
=====

Developed and tested on R 3.0.3 with Windows 7, Windows 8, Ubuntu 13.04 Linux, Mac OS X 10.9.3.

####Setup:
##### Windows
      -[Install R 3.0.3](http://cran.r-project.org/bin/windows/base/old/3.0.3/ "Windows")
##### Mac
      -[Install R-3.0.3.pkg](http://cran.r-project.org/bin/macosx/old/ "Mac OS X")
      -[Install GTK](http://r.research.att.com/libs/GTK_2.24.17-X11.pkg "Mac OS X")
      -Install Xcode from the App Store
      -[Install Xquartz](http://xquartz.macosforge.org/)
##### Linux
Run in terminal: `sudo apt-get install r-base r-cran-rgl libgtk2.0-dev`

2. Load package and dependencies from R Console:
```R
source("http://retistruct.r-forge.r-project.org/install.R") ## retistruct
source("http://retistruct.r-forge.r-project.org/install-gui.R") ## retistruct interface
install.packages('devtools') 
devtools::install_github('bcohn12/retina'); library(retina)
```
3. Demo:
```R
demo('retinaplot')
```


Support
=====
In-package documentation
```R
?retinaplot
```
[Post a github issue ticket (15 seconds)](https://github.com/bcohn12/retina/issues/new "Post an issue ticket")

References
=====
Sterratt et al 2013

GPL-2 License



