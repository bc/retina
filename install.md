Install
=====

Developed and tested on R 3.1.1 with Windows 8.1

####1. Setup:
##### Windows
[Install R 3.1.2](http://cran.r-project.org/bin/windows/base/ "Windows")
[Install RTools](http://cran.r-project.org/bin/windows/Rtools/ "Windows")

##### Linux
Run in terminal: `sudo apt-get install r-base r-cran-rgl libgtk2.0-dev`
##### Mac (experimental)
[Install R 3.1.2](http://cran.r-project.org/bin/macosx/ "Mac OS X")  
[Install GTK](http://r.research.att.com/libs/GTK_2.24.17-X11.pkg "Mac OS X")  
Install Xcode from the App Store  
[Install Xquartz](http://xquartz.macosforge.org/)

####2. Load dependencies from the R console:
```R
install.packages(c("retistruct", "geometry"), repos="http://R-Forge.R-project.org", type="source")
retistruct::retistruct() #On Windows, Accept the 'Install GTK+' popup
install.packages('devtools')
devtools::install_github('bcohn12/retina')
```
####3. Load _retina_ and plot the retina of a reef fish (*Pseudodax moluccanus*):
```R
library(retina)
retinaplot(Pmol_753)
```
####Full walkthrough tutorial to use your own retina
[Link to Tutorial on Github](tutorial.md "Tutorial.md")
