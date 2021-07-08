#!/bin/bash
#Do all the steps
#Note that you may want to run this with nohup
#nohup /bin/bash ./run001.sh >logs/run001.txt 2>&1 &

#Exit on failure
set -e

#Environment variables
. set_env.sh

#Create input files
python scripts/mk_input_files_run001.py

#Mesh generation
python $SRCLOC/geom_mk_msh.py params/mesh/run001.yaml
python $SRCLOC/geom_mk_xml.py params/mesh/run001.yaml
python $SRCLOC/geom_mk_hdf5.py params/mesh/run001.yaml

#Simulations
python2 $SRCLOC/simulator_run.py params/model/run001.yaml

#Post-processing
##TODO