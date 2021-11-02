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

# Install dependencies for R packages
sudo apt-get install libcurl4-openssl-dev libxml2-dev

# Make libraries writeable for active user
sudo chown -R $USER /usr/lib/R/library

# Install DESeq2 R package
Rscript install/deseq2.R
