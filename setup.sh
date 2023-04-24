#!/bin/bash

docker build -t tepic assets/tepic
docker build -t tfprio-python assets/docker/tfprio-python
docker build -t tfprio-r assets/docker/tfprio-r
docker build -t tfprio-angular assets/docker/tfprio-angular