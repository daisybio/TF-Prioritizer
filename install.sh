#!/bin/sh

# Prerequisits
# sudo usermod -a -G staff $USER

if R --version; then
    echo "R already installed."
else
    echo "Installing R."

    # Install necessary dependencies for adding a repo over HTTPS
    sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common

    # Add repo key and repo itself
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

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
sudo apt-get install libcurl4-openssl-dev libxml2-dev

# Make libraries writeable for active user
sudo chown -R $USER /usr/lib/R/library

# Install DESeq2 R package
Rscript install/deseq2.R

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

    echo "PATH=\"\$PATH:$HOME/.meme/libexec/meme-5.4.1:/home/nico/.meme/bin\"" >>$HOME/.profile
    rm -rf meme-5.4.1

    echo "Finished MEME installation."
fi

# Install command line IGV and igvtools
wget https://data.broadinstitute.org/igv/projects/downloads/2.11/IGV_2.11.2.zip
unzip IGV_2.11.2.zip
rm IGV_2.11.2.zip

# Install GCC 9.3.0
sudo apt-get install build-essential