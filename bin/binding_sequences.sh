#!/bin/bash

grep -v "^#" $1 \
    | tail -n +2 \
    | awk -v OFS='\t' '{sub(/^>/,"",$7); gsub(/\([^()]*\)/,"",$1); split($7,loc,/:|-/); print loc[1],loc[2]+$4,loc[2]+$4+$5,$1,$2,$6 >> $1".bed"}'
