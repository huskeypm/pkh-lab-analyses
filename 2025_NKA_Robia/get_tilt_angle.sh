#### This script computes the tilt angle of TM3 relative to the bilayer plane ####

#!/bin/bash

path=/home/xfang2/backup/faust_backup/nka/md_analysis

for sys in wt g303a g303r
do

# traj with only solutes
cat > tmp.trajin << EOD
parm $path/${sys}_noWAT.top parm [top]
trajin $path/${sys}_noWAT.dcd 1 -1 parm [top]
vector VA :257,258,259@N,CA,C,O :280,281,282@N,CA,C,O
vector VM corrplane :PA|:PC|:PE|:PS|:OL|:CHL
vectormath vec1 VM vec2 VA out angle_${sys}.dat name AngleA dotangle
EOD
cpptraj -i tmp.trajin

done

rm -rf tmp.trajin
