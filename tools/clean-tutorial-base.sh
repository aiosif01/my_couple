#!/bin/sh
set -e -u
# -e : non-zero exit status (indicating an error)
# -u : exit if there’s an attempt to use an undefined variable

# shellcheck disable=SC1091
. ../tools/cleaning-tools.sh

clean_tutorial .
clean_precice_logs .
rm -fv ./*.log
rm -fv ./*.vtu

