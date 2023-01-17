#!/bin/sh

./mvnw clean compile assembly:single
mkdir -p bin
cp target/TFPRIO-1.0-bin.jar bin/TF-Prioritizer.jar
docker build . -t ghcr.io/biomedbigdata/tfprio
docker push ghcr.io/biomedbigdata/tfprio