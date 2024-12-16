#!/bin/bash

# Check if a file is provided as an argument
if [ $# -lt 1 ]; then
    echo "Error: Please specify a SAM file as an argument."
    exit 1
fi

# Get the file passed as an argument
file="$1"

# Verify if the argument is a file not a directory
if [ -d "$file" ]; then
    echo "Error: $file is a directory. Please provide a SAM file."
    exit 1
fi

# Check if the file exists and is accessible
if [ ! -f "$file" ]; then
    echo "Error: The file $file does not exist or is not accessible."
    exit 1
fi

# Check if the file contains the required SAM headers starting with '@'
if ! grep -q "^@" "$file"; then
    echo "Error: The file $file does not contain valid SAM headers."
    exit 1
fi

# If all checks pass, execute the Python script
python3 Systeme.py "$file"

# Display a success message if everything is validated
if [ $? -eq 0 ]; then
    echo "Veracity test passed: The file $file is valid."
else
    echo "Error: An issue occurred while executing the Python script."
fi

