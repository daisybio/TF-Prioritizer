#!/bin/bash

docker build -t tepic assets/tepic
docker build -t tfprio-python assets/docker/tfprio-python
docker build -t dynamite assets/docker/dynamite