#!/bin/bash

#calculate the CaM/CaMBR distance, CaM's radius of gyration and RMSD to native complex

cwd=`pwd`

for z in 0M 0.15M 0.5M
do
for i in `seq 1 200`
do
cd $cwd/$z/$i

cat > trajin << EOD
trajin viz.pdb
reference 4q5u_CG.pdb [reff]
rms test1 ref [reff]  out rmsd.dat
rms test2 ref [reff] @1-140 out Nrmsd.dat
rms test3 ref [reff] @160-300 out Crmsd.dat
rms test4 ref [reff] @1-300 out CaMrmsd.dat
rms test5 ref [reff] @301-354 out CaMBRrmsd.dat
radgyr @1-300 out RoG.dat mass nomax
distance d1 @1-300 @301-354 out dis.dat
EOD

cp $cwd/4q5u_CG.pdb .
cp $cwd/trajin .
cp $cwd/top.pdb .

cpptraj -p top.pdb -i trajin
done
done

