#!/bin/bash

# Script to check ORCA output files and move failed jobs to a failed directory

# Find all .out files in the current directory and subdirectories
find . -name "*.out" -type f | while read -r outfile; do
    # Get the last 5 lines of the file
    last_lines=$(tail -n 5 "$outfile")
    
    # Check if "****ORCA TERMINATED NORMALLY****" is in the last 5 lines
    if ! echo "$last_lines" | grep -q "ORCA TERMINATED NORMALLY"; then
        # Get the directory containing the .out file
        dir=$(dirname "$outfile")
        
        # Get the base filename without extension
        base=$(basename "$outfile" .out)
        
        # Create failed directory in the same location as the .out file
        failed_dir="$dir/failed"
        mkdir -p "$failed_dir"
        
        echo "Moving failed job: $base from $dir"
        
        # Move all files matching the pattern {base}.* to the failed directory
        find "$dir" -maxdepth 1 -name "${base}.*" -type f -exec mv {} "$failed_dir/" \;
    fi
done

echo "Done checking ORCA jobs."
