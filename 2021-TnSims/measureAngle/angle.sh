#!/bin/bash


# we use the VGM program to measure the angle between the
# two helics of the EF-hands

# Also since the VGM program only reads the PDB format, we will need to use 
# CPPTRAJ to extract some frames from the trajectory and save them as PDB files

# Based on the manuscript, the ref pdb is residue 14-25 and 38-47 from 1ap4.pdb



# the program cpptraj and vgm are preinstalled in our workstation

cp header anginfo

cat > trajin << EOD
trajin r1prod1.nc 1 -1 20
trajin r1prod2.nc 1 -1 20
trajin r1prod3.nc 1 -1 20
trajin r2prod1.nc 1 -1 20
trajin r2prod2.nc 1 -1 20
trajin r2prod3.nc 1 -1 20
trajin r2prod4.nc 1 -1 20
trajin r3prod1.nc 1 -1 20
trajin r3prod2.nc 1 -1 20
strip !(:14-25,38-47)
trajout ef.pdb
EOD

cpptraj -p struc.pdb -i trajin  

nmoles=`grep MODEL ef.pdb  | wc -l| awk '{ print $1 }'`
lines=365 

grep -n MODEL ef.pdb  | awk -F ":" '{ print $1 }' > list

for i in `seq 1 $nmoles`
do
	
	startline=`awk '{ if ( NR == '$i' ) print }' list`
	endline=$((startline + lines))
	awk '{ if ( NR > '$startline' && NR <= '$endline' ) print }' ef.pdb > temp.pdb
	vgm input_file | tail -n 1 >> anginfo
done




