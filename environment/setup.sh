#!/bin/bash

export DEBIAN_FRONTEND=noninteractive

DIRECTORY=$(dirname "$0")

apt-get install -y \
  openjdk-17-jre-headless \
  r-base \
  python3 \
  python3-pip \
  cmake \
  bedtools \
  unzip \
  xvfb

Rscript "$DIRECTORY"/setup.R
python3 -m pip install -r "$DIRECTORY"/python_dependencies.txt
