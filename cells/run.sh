#!/bin/bash
set -e -u

. ../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# Temporarily disable unbound variable errors
set +u
source /home/ale/biodynamo-v1.05.124/bin/thisbdm.sh
set -u

bdm run

close_log
