#!/bin/bash

# Navigate into the data directory
cd data || exit 1

# Loop through all matching barcode files
for file in Pla_*_barcodes.tsv.gz; do
    # Extract sample ID (everything up to the second underscore)
    sample_id=$(echo "$file" | cut -d'_' -f1-2)

    # Create a folder for the sample
    mkdir -p "$sample_id"

    # Move and rename associated files
    for suffix in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
        src_file="${sample_id}_${suffix}"
        dst_file="$sample_id/$suffix"

        if [ -f "$src_file" ]; then
            mv "$src_file" "$dst_file"
        fi
    done
done

