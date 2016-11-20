![retina_logo](https://cloud.githubusercontent.com/assets/4623063/8342959/8206dd04-1a85-11e5-8d00-d58866c99d66.jpg)
======
Open-source tools for visualizing and comparing retinal cell density data.

#Project description

Biologists who study vision are fascinated by the retina, the rear surface of the eye which collects light and transmits an 'image' to the brain. Through evolution, eyes have been exceptionally diverse across the animal kingdom, yielding a wide variety of shapes, sizes, colors, and chemical processes. For this reason, the eye is often used to understand how complex traits evolve with ecology, the interactions of the animal with its environment. A retinal ganglion cell map helps biologists visualize receptive cells across the rear surface of the vertebrate eye. This topographic map shows density in a way which highlights areas of higher visual performance, and this can be used to understand what regions of an animal's visual field are most clear-in-view. For a human, we only have one point on the eye where things are clear- that's why we can really only see a few words at a time while reading, for example. Some fish have multiple centers of high visual performance, more formally referred to as a peak of retinal acuity.

This project provides an accurate and rigorous tool for visualizing data from an animal's retina, the rear surface of the eye. In particular, our open-source package can be quickly installed by vision researchers worldwide using completely free software, to begin generating meaningful comparisons between different eyes. In the past, the most frequently used way to visualize data on the retina was to slice and flatten the retina onto a microscope slide, then present the data in that form. These _retinal maps_ can take several hours to construct by hand or with proprietary software, and have been the industry standard for over thirty five years, in hundreds of publications. With our software, a researcher can produce many maps in one hour, make consistent comparisons, and generate useful statistics. The software is also customizable to meet the user's specific requirements, as it is open-source. We hope our software will be useful for vision scientists and hope to receive feedback that results in further additions to the program.

Linux: [![Travis-CI Build Status](https://travis-ci.org/bcohn12/retina.png?branch=master)](https://travis-ci.org/bcohn12/retina)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/v7vav80absnsh9jf?svg=true)](https://ci.appveyor.com/project/bcohn12/retina)

<img src="tutorial_pix/retina_plot_output_pmol753.jpg" width=400 alt="some_text">

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


Development Notes
=====
run R from the retina directory after cloning
```
library(devtools)
```

Press
=====
[Keck Science Department Announcement](http://www.kecksci.claremont.edu/News/Newsdetail.asp?NewsID=92 "KSD")

