#!/bin/bash

# Prompt the user for input of the specific main.c file to use
read -p "Enter the specific main.c file to use (e.g., 4-Feature-WT-ACK-Main.c): " main_file

# Prompt the user for input of the output binary name
read -p "Enter the name of the output binary (e.g., 4-Feat-WT-ACK): " output

# Define the directory where the .c files are located
source_dir="/Users/pwoods/week15August2024/4-Feature-Testing/4FeatureCodeTest"  # Replace with the actual directory

# Find all .c files in the source directory, excluding all files ending with *Main.c
other_files=$(find "$source_dir" -type f -name "*.c" ! -name "*Main.c")

# Compile the specified main.c file with other .c files
gcc -o "$output" "$source_dir/$main_file" $other_files

# Inform the user
echo "Compiled $main_file with other source files into $output."

