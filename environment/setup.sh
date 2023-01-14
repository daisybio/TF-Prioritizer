#!/bin/bash

export DEBIAN_FRONTEND=noninteractive

DIRECTORY=$(dirname "$0")

apt-get update && apt-get install -y \
  openjdk-17-jre-headless \
  r-base \
  python3 \
  python3-pip \
  cmake \
  bedtools \
  unzip \
  xvfb \
  wget \
  libcurl4-openssl-dev libxml2-dev r-cran-httr

mkdir -p "$RGTDATA"
python3 -m pip install -r "$DIRECTORY"/python_dependencies.txt
chmod -R 777 "$RGTDATA"

Rscript "$DIRECTORY"/setup.R
