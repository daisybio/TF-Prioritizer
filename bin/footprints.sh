#!/bin/bash

input="$1"
output="$2"
search_type="$3"
max_dist="$4"

if [[ $search_type == "incl_between" ]]; then
    start_col=2
    end_col=6
elif [[ $search_type == "excl_between" ]]; then
    start_col=3
    end_col=5
else 
    printf "Invalid search type\n" >&2
    exit 1
fi

# Simplify and sort input
cut -f1-3 "$input" | bedtools sort > "sorted.bed"

# 1. Find the closest downstream interval for each interval
# -D a: Add the distance to the closest interval to the end of each line
# -io: Ignore overlapping
# -iu: Ignore upstream

# 2. Filter out intervals that are too far away or where no match was found
# 3. Creat footprints from paired intervals
# 4. Merge footprints (only necessary for incl_between)

bedtools closest -D a -a sorted.bed -b sorted.bed -io -iu \
  | awk -v cutoff=$max_dist '$NF <= cutoff && $NF >= 0' \
  | cut -f1,"$start_col","$end_col" \
  | bedtools merge > "$output"

