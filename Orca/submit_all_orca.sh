#!/bin/bash

# Script to submit ORCA .inp files if corresponding .out file is missing or incomplete

for inp_file in *.inp; do
    # Check if any .inp files exist
    if [ ! -e "$inp_file" ]; then
        echo "No .inp files found in current directory"
        exit 0
    fi
    
    # Get the base name without extension
    base_name="${inp_file%.inp}"
    out_file="${base_name}.out"
    
    submit=false
    
    # Check if .out file doesn't exist
    if [ ! -f "$out_file" ]; then
        echo "Output file $out_file not found. Submitting $inp_file"
        submit=true
    else
        # Check if ORCA terminated normally in the last 5 lines
        if ! tail -n 5 "$out_file" | grep -q "ORCA TERMINATED NORMALLY"; then
            echo "Output file $out_file incomplete. Resubmitting $inp_file"
            submit=true
        else
            echo "Output file $out_file already complete. Skipping $inp_file"
        fi
    fi
    
    # Submit the job if needed
    if [ "$submit" = true ]; then
        suborca -p 2 -m 8000 -v 601 "$inp_file"
    fi
done