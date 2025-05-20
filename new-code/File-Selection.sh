#!/bin/bash

## ---- Initial Verification ----

# Store all .c and .txt files recursively
c_files=$(find . -name "*.c")
txt_files=$(find . -name "*.txt")

# Check for .c files
if [ -z "$c_files" ]; then
  echo "No .c files found."
  exit 1
fi

# Check for .txt files
if [ -z "$txt_files" ]; then
  echo "No .txt files found."
  exit 1
fi

## ---- Initial Prompts and Variables ----
read -p "Project Name: " project_name

# Select Main.c file
echo "Select a Main.c file:"
select file in $c_files; do
  if [ -n "$file" ]; then
    echo "You selected: $file"
    main_file="$file"
    break
  else
    echo "Invalid selection. Please try again."
  fi
done

# Select input .txt file
echo "Select an input (.txt) file:"
select file in $txt_files; do
  if [ -n "$file" ]; then
    echo "You selected: $file"
    input_file="$file"
    break
  else
    echo "Invalid selection. Please try again."
  fi
done

# Choose output directory and binary path
read -p "Enter the desired output directory: " output_dir
read -p "Enter the comparison used, with the formatting 'WT_ACK': " comparison_name
read -p "Enter the path and name for the output binary (e.g. ./lda_analysis): " binary_path
mkdir -p "$output_dir"

# Select feature analysis mode
echo "Select feature analysis mode:"
echo "1. Analyze 1 feature"
echo "2. Analyze 2 features"
echo "3. Analyze 3 features"
echo "4. Analyze 4 features"
echo "5. Analyze features 1-3"
echo "6. Analyze all (1-4)"
read -p "Enter option (1-6): " feature_option

# Create necessary output files
echo "Preparing output files..."

create_output_file() {
  local n_feat="$1"
  local output_file="${output_dir}/${comparison_name}_${n_feat}feature_output.txt"
  touch "$output_file"
  if [[ $? -eq 0 ]]; then
    echo "Created: $output_file"
  else
    echo "Failed to create: $output_file"
    exit 1
  fi
}

case "$feature_option" in
  1) create_output_file 1 ;;
  2) create_output_file 2 ;;
  3) create_output_file 3 ;;
  4) create_output_file 4 ;;
  5)
    for i in {1..3}; do
      create_output_file "$i"
    done
    ;;
  6)
    for i in {1..4}; do
      create_output_file "$i"
    done
    ;;
  *)
    echo "Invalid option selected for output preparation."
    exit 1
    ;;
esac

# ---- Final Verification Before Compilation ----
echo ""
echo "Please verify the following selections:"
echo "Project Name: $project_name"
echo "Main File: $main_file"
echo "Input File: $input_file"
echo "Output Directory: $output_dir"
echo "Binary Path: $binary_path"
echo "Feature Option: $feature_option"
echo "Output Name: $output_file"
echo ""
echo "Would you like to proceed?"
echo "1. Yes"
echo "2. No"
echo "3. Exit"
read -p "> " verify_1

if [[ "$verify_1" = "3" ]]; then
  echo "Exiting."
  exit 1
elif [[ "$verify_1" != "1" ]]; then
  echo "Cancelled. Rerun script to restart."
  exit 1
fi

# ---- Compilation ----
source_dir=$(dirname "$main_file")
other_files=$(find "$source_dir" -type f -name "*.c" ! -name "*Main.c")

gcc -o "$binary_path" "$main_file" $other_files -lm
if [[ $? -ne 0 ]]; then
  echo "Compilation failed. Exiting."
  exit 1
else
  echo "Compilation successful. Binary created: $binary_path"
fi

# ---- Run Binary Prompt ----
echo "Run the binary now?"
echo "1. Yes"
echo "2. No"
read -p "> " run_now

if [[ "$run_now" != "1" ]]; then
  echo "Execution skipped. You may run the binary manually later."
  exit 0
fi

# ---- Run Binary Based on Selected Option ----
run_lda () {
  local input="$1"
  local out_dir="$2"
  local feat="$3"
  local output_file="${out_dir}/${comparison_name}_${feat}feature_output.txt"
  "$binary_path" "$input" "$output_file" "$feat"
  echo "Output saved to: $output_file"

}

case "$feature_option" in
  1) run_lda "$input_file" "$output_dir" 1 ;;
  2) run_lda "$input_file" "$output_dir" 2 ;;
  3) run_lda "$input_file" "$output_dir" 3 ;;
  4) run_lda "$input_file" "$output_dir" 4 ;;
  5)
    for i in {1..3}; do
      run_lda "$input_file" "$output_dir" "$i"
    done
    ;;
  6)
    for i in {1..4}; do
      run_lda "$input_file" "$output_dir" "$i"
    done
    ;;
  *)
    echo "Invalid option during execution."
    exit 1
    ;;
esac

echo "Feature analysis complete. Results saved in: $output_dir"

