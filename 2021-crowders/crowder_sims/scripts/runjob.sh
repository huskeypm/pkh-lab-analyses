#!/bin/bash
#Prepare the job files for the specified analysis, and run it

#Exit on failure
set -e

##source ./set_env.sh
python3 -m simproc requests/mk_${1}.yaml
python3 -m simproc requests/generated/${1}.yaml
