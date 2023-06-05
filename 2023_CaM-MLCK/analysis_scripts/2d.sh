#!/bin/bash

path=/home/xfang2/backup/faust_backup/mlck/smMLCK/md/2_mlck_cam/analyses/shortMD
systems=(cam_t34k cam_dN mcam cam_k30e cam_g40d cam_m36i scam1 scam4)

for sys in ${systems[@]}
do

for comp in 2
do

for rep in {1..50}
do

cat > tmp.trajin << EOD
readdata $path/dist/cam_mlck/dist_${sys}_${comp}_${rep}.dat name set1
readdata $path/dist/dfg_hrd/dist_${sys}_${comp}_${rep}.dat name set2
hist set1 set2 bins 50 out cm_dfg_${sys}_${comp}_${rep}.gnu name PC12 free 310
hist set1 set2 bins 50 out cm_dfg_prob_${sys}_${comp}_${rep}.gnu name PC12p norm
EOD

cpptraj -i tmp.trajin

done
done
done
