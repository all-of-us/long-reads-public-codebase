#!/bin/bash

# Run CAVIAR (default settings) in every per-gene folder written by caviar_inputs.py.
#
# Input:  caviar/<gene>/variant.ld and variant.zscore
# Output: CAVIAR results in each folder, with output prefix "test"
# Requires the CAVIAR binary on PATH as "CAVIAR". README.md section 3, step 4.
#
# Usage: bash run_caviar.sh

base_dir="caviar"

for dir in "$base_dir"/*; do
    if [ -d "$dir" ]; then  
        pushd "$dir" > /dev/null  
        if [[ -f "variant.ld" && -f "variant.zscore" ]]; then
            CAVIAR -l variant.ld -z variant.zscore -o test
            echo "CAVIAR ran successfully in $dir"
        else
            echo "Required files not found in $dir"
        fi
        popd > /dev/null  
    fi
done
