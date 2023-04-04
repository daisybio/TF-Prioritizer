#!/bin/bash

DIRECTORY=$(dirname "$0")

# Prepare node installation
NPM_CACHE=/srv/dependencies/npm
mkdir -p $NPM_CACHE
apt-get update && apt-get install -y curl
curl -sL https://deb.nodesource.com/setup_18.x | bash -

# bedtools: Required by TEPIC
# unzip: Required for ChipAtlas download extraction
# xvfb: Required for headless IGV
# wget: Required for HINT genome download
# curl: Required for angular installation
# libcurl4-openssl-dev, libxml2-dev, r-cran-httr: Required by R packages
apt-get update && apt-get install -y \
  nodejs \
  openjdk-17-jre-headless \
  openjdk-17-jdk \
  r-base \
  python3 \
  python3-pip \
  cmake \
  bedtools \
  samtools \
  unzip \
  xvfb \
  wget \
  curl \
  pngquant \
  libcurl4-openssl-dev libxml2-dev r-cran-httr

npm install npm@latest -g
npm install -g @angular/cli
npm config set cache $NPM_CACHE --global
chmod -R 777 $NPM_CACHE

# Setup matplotlib cache directory
mkdir -p "$MPLCONFIGDIR"
chmod -R 777 "$MPLCONFIGDIR"

# Install python packages, RGTDATA=Cache directory for HINT
mkdir -p "$RGTDATA"
python3 -m pip install -r "$DIRECTORY"/python_dependencies.txt
chmod -R 777 "$RGTDATA"

# Install R packages
Rscript "$DIRECTORY"/setup.R

# Install IGV
mkdir -p "$IGV"
wget https://data.broadinstitute.org/igv/projects/downloads/2.15/IGV_2.15.4.zip -O "$IGV"/IGV.zip
unzip "$IGV"/IGV.zip -d "$IGV"
mv "$IGV"/IGV_2.15.4/* "$IGV"
rm -rf "$IGV"/IGV_2.15.4 "$IGV"/IGV.zip

mkdir -p "$IGV_CACHE"
echo "PORT_ENABLED=false" > "$IGV_CACHE"/prefs.properties
chmod -R 777 "$IGV_CACHE"