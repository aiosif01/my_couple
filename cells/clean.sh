#!/bin/bash

# Script to clean up BioDynaMo generated files in the cells directory

echo "Cleaning BioDynaMo simulation output..."

# Remove log files
echo "Removing log files..."
rm -f cells.log

# Clean output directory but keep the directory itself
if [ -d "output" ]; then
  echo "Cleaning output directory..."
  rm -rf output/*
fi

# Remove precice-profiling directory contents
if [ -d "precice-profiling" ]; then
  echo "Cleaning precice-profiling directory..."
  rm -rf precice-profiling/*
fi

# Clean build artifacts if needed
if [ -d "build" ]; then
  echo "Cleaning build directory..."
  # Only remove content inside build folder, not the directory itself
  rm -rf build/*
fi

# Remove any temporary files
echo "Removing temporary files..."
rm -f *~
rm -f *.dat
rm -f core.*

echo "Clean up complete."