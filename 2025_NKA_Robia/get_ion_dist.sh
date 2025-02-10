#### This script computes distance of the Na+ ions from the binding site as a quantification of ion release ####


#!/bin/bash

path=/home/xfang2/backup/faust_backup/nka/md_analysis

for sys in wt g303a g303r wt_dimer g303a_dimer g303r_dimer
do

# traj with only solutes
cat > tmp.trajin << EOD
parm $path/${sys}.pdb parm [top]
trajin $path/${sys}.dcd 1 -1 parm [top]
EOD

if [ $sys == wt_dimer ] || [ $sys == g303a_dimer ] || [ $sys == g303r_dimer ]
then
cat >> tmp.trajin << EOD
distance dist1 :2669 :302@O|:304@O|:306@OE1|:783@OD1|:787@OD1 out dist_${sys}_Na1_dimer_1.dat
distance dist2 :2670 :302@O|:755@OD1|:751@O|:787@OD1,OD2 out dist_${sys}_Na2_dimer_1.dat
distance dist3 :2671 :902@OE1|:750@O|:905@OD1,OD2 out dist_${sys}_Na3_dimer_1.dat
distance dist4 :2679 :1636@O|:1638@O|:1640@OE1|:2117@OD1|:2121@OD1 out dist_${sys}_Na1_dimer_2.dat
distance dist5 :2680 :1636@O|:2089@OD1|:2085@O|:2121@OD1,OD2 out dist_${sys}_Na2_dimer_2.dat
distance dist6 :2681 :2236@OE1|:2084@O|:2239@OD1,OD2 out dist_${sys}_Na3_dimer_2.dat
EOD

else
cat >> tmp.trajin << EOD
distance dist1 :1335 :302@O|:304@O|:306@OE1|:783@OD1|:787@OD1 out dist_${sys}_Na1.dat
distance dist2 :1336 :302@O|:755@OD1|:751@O|:787@OD1,OD2 out dist_${sys}_Na2.dat
distance dist3 :1337 :902@OE1|:750@O|:905@OD1,OD2 out dist_${sys}_Na3.dat
EOD

fi

cpptraj -i tmp.trajin

done

rm -rf tmp.trajin
