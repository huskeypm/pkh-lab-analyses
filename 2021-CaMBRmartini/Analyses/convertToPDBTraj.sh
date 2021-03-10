#!/bin/bash

# convert gromacs .xtc traj to pdb format with bonds added
# for visualization purpose

cwd=`pwd`
for z in 0M 0.15M 0.5M
do
for i in `seq 1 200`
do
cd $cwd/$z/$i

if [ -f viz.pdb ]; then
        continue
else
        echo $i

# generate a pdb to serve as the topology for cpptraj
gmx editconf -f equi.gro -o struc.pdb

# cpptraj to autoimage
cat > trajin << EOD
trajin prod.xtc 1 -1 10
autoimage
trajout newprod.xtc
EOD

cpptraj -p struc.pdb -i trajin 2> aa.log

echo 1 1 | gmx trjconv -f newprod.xtc -s prod.tpr -fit rot+trans -o viz.pdb -conect 2> bblog
fi
done
done

