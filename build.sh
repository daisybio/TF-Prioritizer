#!/bin/sh

./mvnw clean compile assembly:single
cp target/TFPRIO-1.0-jar-with-dependencies.jar bin/TF-Prioritizer.jar
docker build . -t ghcr.io/biomedbigdata/tfprio
docker push ghcr.io/biomedbigdata/tfprio