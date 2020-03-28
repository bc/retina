git clone git@github.com:bc/retina.git -o StrictHostKeyChecking=no
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install r-base libssl-dev -y
sudo apt-get install libcurl4-openssl-dev libxml2-dev

