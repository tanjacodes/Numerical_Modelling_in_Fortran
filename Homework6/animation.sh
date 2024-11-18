#!/bin/bash

# Define the source and destination directories
SOURCE_DIR="/home/tanja/Documents/Studium/Master/'1. Semester'/Numerical_Modelling_in_Fortran/Homework6"
DEST_DIR="/mnt/c/Users/41766/Documents/1RW/5/Numerical Modelling in Fortran/Homework/Week6/Data"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Move the .dat files from the source to the destination
for i in $(seq -w 0 200); do
    mv "$SOURCE_DIR/T_output_0$i.dat" "$DEST_DIR/"
    mv "$SOURCE_DIR/S_output_0$i.dat" "$DEST_DIR/"
done

echo "Files moved successfully from $SOURCE_DIR to $DEST_DIR."
