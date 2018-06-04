# Install

Developed and tested on R 3.1.1 with Windows 8.1

## 1\. Setup:

### Windows

[Install R 3.1.2](http://cran.r-project.org/bin/windows/base/ "Windows") [Install RTools](http://cran.r-project.org/bin/windows/Rtools/ "Windows")

### Linux

Run in terminal:

```bash
sudo apt-get -y build-dep libcurl4-gnutls-dev
sudo apt-get -y install libcurl4-gnutls-dev
sudo apt-get -y install libglu1-mesa-dev
sudo apt-get -y install r-base-core r-base r-cran-rgl libgtk2.0-dev
```

### Mac

[Install R 3.5.0](http://cran.r-project.org/bin/macosx/ "Mac OS X")<br>
Run in the terminal:

```bash
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install gtk+
export PKG_CONFIG_PATH=/usr/X11/lib/pkgconfig:$PKG_CONFIG_PATH
brew install caskroom/cask/brew-cask
brew cask install xquartz
```

Run in R

```r
install.packages(c("cairoDevice", "RGtk2"), type="source")
install.packages('devtools')
devtools::install_github('bc/retina')
```

Finally, install Xcode from the App Store

## 2\. Load dependencies from the R console:

```r
install.packages('retistruct') #if you have an issue, comment on it here to get it fixed: https://github.com/davidcsterratt/retistruct/issues/new
retistruct::retistruct() #On Windows, Accept the 'Install GTK+' popup
install.packages('devtools')
devtools::install_github('bc/retina')
```

## 3\. Load _retina_ and plot the retina of a reef fish (_Pseudodax moluccanus_):

```r
library(retina)
sample_retina <- banded_gecko()
retinaplot(sample_retina)
```

## Full walkthrough tutorial to use your own retina

[Link to Tutorial on Github](tutorial.md "Tutorial.md")
