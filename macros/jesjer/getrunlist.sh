#!/bin/bash

# Check for input file
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <file_with_paths.txt>"
    exit 1
fi

input_file="$1"
run_numbers=()

while IFS= read -r line; do
    # Extract the last part of the filename
    filename=$(basename "$line")
    
    # Extract the run number using pattern matching
    if [[ "$filename" =~ -([0-9]{8})\.root$ ]]; then
        full_number="${BASH_REMATCH[1]}"
        # Take the last 5 digits of the 8-digit number
        run_number="${full_number: -5}"
        run_numbers+=("$run_number")
    fi
done < "$input_file"

# Output as comma-separated list
IFS=','; echo "${run_numbers[*]}"
