#!/bin/bash


# Converts a bed file into a gff file
# $1 is the input bed file
# $2 is the output gff file

awk 'BEGIN {OFS="\t"} {print $1, "bed2gff", "feature", $2+1, $3, ".", $6, ".", "ID=" $4}' $1 > $2
