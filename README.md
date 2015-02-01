Linux:[![Travis-CI Build Status](https://travis-ci.org/bcohn12/retina.png?branch=master)](https://travis-ci.org/bcohn12/retina), Windows:[![Build status](https://ci.appveyor.com/api/projects/status/v7vav80absnsh9jf?svg=true)](https://ci.appveyor.com/project/bcohn12/retina)


retina in R
======

Open-source tools for visualizing and comparing retinal cell density data.
<img src="tutorial_pix/retina_plot_output_pmol753.jpg" width=400 alt="some_text">

> A retinal ganglion cell map helps biologists visualize receptive cells across the rear surface of the vertebrate eye. This topographic map shows density in a way which highlights areas of higher visual performance.

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

[Link to installation guide](install.md "Installation Page")

####Load _retina_ and  Pseudodax Moluccanus:
```R
library(retina) #see the installation guide first
retinaplot(Pmol_753)
```
####Full walkthrough tutorial to use your own retina
[Link to Tutorial on Github](tutorial.md "Tutorial.md")

Support
=====
In-package documentation
```R
?retinaplot
```
[Post a github issue ticket (15 seconds)](https://github.com/bcohn12/retina/issues/new "Post an issue ticket")

References
=====
1. Sterratt et al 2013

GPL-2 License



