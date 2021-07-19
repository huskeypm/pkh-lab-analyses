#!/bin/bash
#Do all the steps
#Note that you may want to run this with nohup
#nohup /bin/bash ./run004.sh >logs/run004.txt 2>&1 &

#Info
basename="run005"
echo "dir:" `pwd`
echo "node:" `uname -n`
echo "stack size:" `ulimit -s`
echo "w:"
w
echo "date:" `date`
echo "basename:" $basename
echo "---===---"

#Exit on failure
set -e

#Environment variables
echo `date` 'Setting environment variables'
. set_env.sh

#Create input files
echo `date` 'Creating input files'
python scripts/mk_input_files_$basename.py

#Mesh generation
echo `date` 'Mesh generation'
python $SRCLOC/geom_mk_msh.py params/mesh/$basename.yaml
python $SRCLOC/geom_mk_xml.py params/mesh/$basename.yaml
python $SRCLOC/geom_mk_hdf5.py params/mesh/$basename.yaml

#Simulations
echo `date` 'Starting simulations'
simsfile="${basename}_sims.sh"
. $simsfile
#for batch jobs, would need script to check for completion before closing
echo `date` 'Waiting for job completion.'
wait
echo `date` 'Jobs complete.'

#Post-processing
##TODO

#Info
echo "---===---"
echo "date: " `date`
echo "uptime:" `uptime`
