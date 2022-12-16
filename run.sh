#!/bin/sh

./build.sh
java -jar target/TFPRIO-1.0-jar-with-dependencies.jar -o /home/nico/Data/TFPRIO/Runs/test -c /home/nico/Software/TF-Prioritizer/configTemplates/template.json -t 5