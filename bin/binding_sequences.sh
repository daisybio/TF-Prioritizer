#!/bin/bash

INPUT=$1
CUTOFF=$2

grep -v "^#" $INPUT \
    | tail -n +2 \
    | awk "\$2 >= $CUTOFF" \
    | awk -v OFS='\t' '{sub(/^>/,"",$7); gsub(/\([^()]*\)/,"",$1); split($7,loc,/:|-/); print loc[1],loc[2]+$4,loc[2]+$4+$5,$2 >> $1".unsorted"}'


# Sort all the bed files
for i in *.unsorted; do bedtools sort -i $i > $(basename $i .unsorted).bed; done

# Delete the unsorted files
rm -f *.unsorted