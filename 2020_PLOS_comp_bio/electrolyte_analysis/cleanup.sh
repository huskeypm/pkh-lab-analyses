#!/bin/bash

#Exit on failure
set -e

#Get basename from command line
if [ ${#1} = 0 ]
then
  echo "Must provide basename"
  exit 1
fi
basename=$1

#Remove log files
rm -f ./logs/${basename}*.txt

#Remove mesh files
rm -Rf ./mesh/geo/$basename
rm -Rf ./mesh/gmsh_out/$basename
rm -Rf ./mesh/metadata/$basename
rm -Rf ./mesh/msh/$basename
rm -Rf ./mesh/xml/$basename
rm -Rf ./mesh/hdf5/$basename

#Remove solution files
rm -Rf ./solutions/$basename

#Remove parameter files
rm -f ./params/mesh/${basename}*.yaml
rm -f ./params/model/${basename}*.yaml

#Remove codebook
rm -f ./codebook_$basename.yaml

#Remove simulation scripts
rm -f ./${basename}_sims.sh
