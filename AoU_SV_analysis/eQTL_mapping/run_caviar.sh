#!/bin/bash

base_dir="caviar"

for dir in "$base_dir"/*; do
    if [ -d "$dir" ]; then  
        pushd "$dir" > /dev/null  
        if [[ -f "variant.ld" && -f "variant.zscore" ]]; then
            CAVIAR -l variant.ld -z variant.zscore -o test -r 0.95 -c 1 -f 1
            echo "CAVIAR ran successfully in $dir"
        else
            echo "Required files not found in $dir"
        fi
        popd > /dev/null  
    fi
done
