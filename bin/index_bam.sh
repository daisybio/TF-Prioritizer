#!/bin/bash

# Indexes the bam file (apply after reheadering and re-sorting)
# $1 is the path to the file
# Output file will be in the same path as the input file with ".bai" attached

samtools index "$1"
