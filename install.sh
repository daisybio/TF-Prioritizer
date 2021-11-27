#!/bin/sh

if R --version; then
    echo "R already installed."
else
    echo "Installing R."

    # Install necessary dependencies for adding a repo over HTTPS
    sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common

    # Add repo key and repo itself
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

    # Install libicu66
    wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu66_66.1-2ubuntu2_amd64.deb
    sudo apt update && sudo dpkg -i libicu66_66.1-2ubuntu2_amd64.deb
    rm libicu66_66.1-2ubuntu2_amd64.deb

    # Install R
    sudo apt-get install r-base

    sudo apt-get install build-essential

    echo "Finished R installation."
fi

if java --version; then
    echo "Java already installed."
else
    echo "Installing Java."

    sudo apt-get install openjdk-17-jdk openjdk-17-demo openjdk-17-doc openjdk-17-jre-headless openjdk-17-source openjdk-17-jre

    echo "Finished Java installation."
fi

# Install dependencies for R packages
sudo apt-get install libcurl4-openssl-dev libxml2-dev libssl-dev

# Make libraries writeable for active user
sudo chown -R "$USER" /usr/local/lib/R/site-library

# Install DESeq2 R package
Rscript install/r_dependencies.R

# Install R packages required for TEPIC
Rscript ext/TEPIC/TEPIC/Code/installRpackages.R

# Install some more linux packages
sudo apt-get install bedtools python3 unzip

# Install required python packages
python3 -m pip install -r install/requirements.txt

# Install MEME Suite
if meme -version; then
    echo "MEME already installed."
else
    echo "Installing MEME."
    wget https://meme-suite.org/meme/meme-software/5.4.1/meme-5.4.1.tar.gz
    tar -xf meme-5.4.1.tar.gz
    rm meme-5.4.1.tar.gz

    cd meme-5.4.1
    ./configure --prefix=$HOME/.meme --enable-build-libxml2 --enable-build-libxslt
    make
    make test
    make install
    cd ..

    touch $HOME/.profile

    echo "PATH=\"\$PATH:$HOME/.meme/libexec/meme-5.4.1:$HOME/.meme/bin\"" >>"$HOME"/.profile
    rm -rf meme-5.4.1

    echo "Finished MEME installation."
fi

# Install command line IGV and igvtools
if [ ! -d "IGV_2.11.2" ]; then
    wget https://data.broadinstitute.org/igv/projects/downloads/2.11/IGV_2.11.2.zip
    unzip IGV_2.11.2.zip
    rm IGV_2.11.2.zip
fi

# Install GCC 9.3.0
sudo apt-get install build-essential

# Install additional Java packages
mkdir -p ext/lib

# Install org.apache.commons.compress
if [ ! -d "lib/commons-compress-1.21" ]; then
    wget https://dlcdn.apache.org//commons/compress/binaries/commons-compress-1.21-bin.tar.gz
    tar -xf commons-compress-1.21-bin.tar.gz
    rm commons-compress-1.21-bin.tar.gz
    mv commons-compress-1.21 ext/lib/commons-compress-1.21
fi

# Install org.apache.commons.cli
if [ ! -d "lib/commons-cli-1.5.0" ]; then
    wget https://dlcdn.apache.org//commons/cli/binaries/commons-cli-1.5.0-bin.tar.gz
    tar -xf commons-cli-1.5.0-bin.tar.gz
    rm commons-cli-1.5.0-bin.tar.gz
    mv commons-cli-1.5.0 ext/lib/commons-cli-1.5.0
fi

# Install weka classifier
if [ ! -d "lib/weka" ]; then
    wget http://www.java2s.com/Code/JarDownload/weka/weka.jar.zip
    unzip weka.jar.zip
    rm weka.jar.zip
    mkdir -p ext/lib/weka
    mv weka.jar ext/lib/weka/weka.jar
fi