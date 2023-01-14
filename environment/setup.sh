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
  libcurl4-openssl-dev libxml2-dev r-cran-httr

Rscript "$DIRECTORY"/setup.R
python3 -m pip install -r "$DIRECTORY"/python_dependencies.txt
