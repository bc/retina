retina Package for R
======

Open-source tools for visualizing and comparing retinal cell density data.


> A retinal ganglion cell map helps biologists visualize the distribution of receptive cells across the rear surface of the vertebrate eye. This topographic map shows density in a way which highlights areas of higher visual performance.


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

1. [Install R from R-Project](http://www.r-project.org/ "R Project Homepage")
The retina package is compatible with 3.1.0 on Windows 7 and Windows 8

2. Load package and dependencies from R Console:
```R
source("http://retistruct.r-forge.r-project.org/install.R") ## retistruct
source("http://retistruct.r-forge.r-project.org/install-gui.R") ## retistruct interface
install.packages('devtools')
devtools::install_github('bcohn12/retina'); library(retina)
```
3. Demo the plots:
```R
demo('retinaplot')
```
required:R 3.1
Retistruct and gui
retina package
run demo

Support and Suggestions
=====
[Post a github issue ticket (15 seconds)](https://github.com/bcohn12/retina/issues/new "Post an issue ticket")


Package documentation (PDF link)

References
=====
1. Sterratt et al 2013

GPL-2 License



