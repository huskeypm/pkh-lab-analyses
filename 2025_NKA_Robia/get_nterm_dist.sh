#### This script computes the N-term distance between two monomers ####

#!/bin/bash

path=/home/xfang2/backup/faust_backup/nka/md_analysis

for sys in wt_dimer g303a_dimer g303r_dimer
do

cat > tmp.trajin << EOD
parm $path/${sys}.pdb parm [top]
trajin $path/${sys}.dcd 1 -1 parm [top]
distance dist1 :996-1000@CA,C,N,O :2330-2334@CA,C,N,O out dist_${sys}.dat
EOD

cpptraj -i tmp.trajin

done


rm -rf tmp.trajin



