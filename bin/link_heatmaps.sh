#!/bin/bash

input="$1"
output="$2"

# Iterate all png files

for heatmap in "$input"/*.png; do
    # Get the name of the heatmap
    name=$(basename "$heatmap" .png)
    
    # Base names are structured like this: Pairing_HM_TF
    # We want to extract Pairing, HM, and TF
    pairing=$(echo "$name" | cut -d'_' -f1)
    hm=$(echo "$name" | cut -d'_' -f2)
    tf=$(echo "$name" | cut -d'_' -f3)

    echo "$pairing"
    echo "$hm"
    echo "$tf"
    echo "Test"

    # Create a new directory: output/TF_heatmaps/pairing
    current_dir="$output/${tf}_heatmaps/$pairing"
    mkdir -p "$current_dir"

    # Softlink the heatmap into the new directory
    ln -s "../../$heatmap" "$current_dir/$name.png"
done