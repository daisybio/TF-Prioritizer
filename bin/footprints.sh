#!/bin/bash

input="$1"
output="$2"
search_type="$3"
max_dist="$4"

cut -f1-3 "$input" | bedtools sort | bedtools merge > preprocessed.bed

# 1. Find the closest downstream interval for each interval
# -D a: Add the distance to the closest interval to the end of each line
# -io: Ignore overlapping
# -iu: Ignore upstream

# 2. Filter out intervals that are too far away or where no match was found
# 3. Creat footprints from paired intervals
# 4. Merge footprints (only necessary for incl_between)

echo "Finding gaps between intervals..."

bedtools closest -D a -a preprocessed.bed -b preprocessed.bed -io -iu \
| awk -v cutoff="$max_dist" '{ if ($NF <= cutoff && $NF >= 0) {print $1 "\t" $3 "\t" $5} }' > temp.bed


if [[ $search_type == "incl_between" ]]; then
    echo "Adding original intervals to footprints..."
    cat temp.bed preprocessed.bed | bedtools sort | bedtools merge > "$output"
elif [[ $search_type == "excl_between" ]]; then
    ln -s temp.bed "$output"
else 
    printf "Invalid search type\n" >&2
    exit 1
fi