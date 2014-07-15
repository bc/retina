Install
=====

Developed and tested on R 3.0.3 with Windows 7, Windows 8, Ubuntu 13.04 Linux, Mac OS X 10.9.3.

####Setup:
##### Windows
[Install R 3.0.3](http://cran.r-project.org/bin/windows/base/old/3.0.3/ "Windows")
##### Mac
[Install R-3.0.3.pkg](http://cran.r-project.org/bin/macosx/old/ "Mac OS X")  
[Install GTK](http://r.research.att.com/libs/GTK_2.24.17-X11.pkg "Mac OS X")  
Install Xcode from the App Store  
[Install Xquartz](http://xquartz.macosforge.org/)
##### Linux
Run in terminal: `sudo apt-get install r-base r-cran-rgl libgtk2.0-dev`

####Open R and load dependencies:
```R
source("http://retistruct.r-forge.r-project.org/install.R") ## retistruct
source("http://retistruct.r-forge.r-project.org/install-gui.R") ## retistruct interface
install.packages('devtools') 
devtools::install_github('bcohn12/retina')
```
####Load _retina_ and  Pseudodax Moluccanus:
```R
library(retina)
retinaplot(Pmol_753)
```
####Full walkthrough tutorial to use your own retina
[Link to Tutorial on Github](tutorial.md "Tutorial.md")
