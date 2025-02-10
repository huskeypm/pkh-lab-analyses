#### This script extracts water coordinates from the MD trajectories ####

#!/bin/bash

path=/home/xfang2/backup/faust_backup/nka/md_analysis

for sys in wt g303a g303r wt_dimer g303a_dimer g303r_dimer
do

# traj with only solutes
cat > tmp.trajin << EOD
parm $path/${sys}.pdb parm [top]
trajin $path/${sys}.dcd 1 -1 30 parm [top]
reference $path/${sys}.pdb
autoimage
EOD

if [ $sys == wt_dimer ] || [ $sys == g303a_dimer ] || [ $sys == g303r_dimer ]
then
cat >> tmp.trajin << EOD
center mass origin :2672-2678|:2692-4185 
align :2672-2678|:2692-4185
strip !(:WAT@O)
trajout coords.pdb pdb
EOD

else
cat >> tmp.trajin << EOD
center mass origin :1339-2517
align :1339-2517
strip !(:WAT@O)
trajout coords.pdb pdb
EOD
fi

cpptraj -i tmp.trajin

# remove useless lines
grep -Ev "TER|MODEL|CRYST|END" coords.pdb > tmp.pdb

# extract only X Y Z corrdinates
awk '{print substr($0, 31, 8), substr($0, 39, 8), substr($0, 47, 8)}' tmp.pdb > ${sys}_coords.dat


done

rm -rf tmp.trajin
