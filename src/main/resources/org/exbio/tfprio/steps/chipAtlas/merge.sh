#!/bin/sh

for f in "$1"/*.bed; do
    TF=$(basename "$f" | cut -d_ -f1)
    awk -v TF="$TF" '{gsub(/^chr/,"",$1);print $1"\t"$2"\t"$3"\t"TF}' "$f"
done | grep -v "^track" | sort -k1,1V -k2,2n | bedtools merge -i stdin -c 4 -o collapse -delim "|" > "$2"
