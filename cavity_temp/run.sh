#!/bin/bash
set -e -u

. ../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# Clean preCICE connection files before running
echo "Cleaning up preCICE connection files..."
rm -rf ../precice-run

foamListTimes -rm
blockMesh
topoSet

../tools/run-openfoam.sh "$@"
. ../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
