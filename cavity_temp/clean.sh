#!/bin/bash

# Script to clean up OpenFOAM generated files in cavity_temp

# Remove all time directories except for 0
echo "Removing time directories..."
foamListTimes -rm

# Remove log files
echo "Removing log files..."
rm -f cavity_temp.log

# Remove precice-profiling directory contents
if [ -d "precice-profiling" ]; then
  echo "Cleaning precice-profiling directory..."
  rm -rf precice-profiling/*
fi

# Remove other common OpenFOAM generated files
echo "Removing other generated files..."
rm -rf processor*/ # For parallel runs
rm -f *.dat *.vtk # Data files
rm -f *.csv # CSV output
rm -f *~  # Temporary files
rm -f log.* # OpenFOAM log files
rm -f ./*.foam # ParaView launcher files

echo "Clean up complete."