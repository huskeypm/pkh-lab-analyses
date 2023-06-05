#!/bin/bash

path=/home/xfang2/backup/faust_backup/mlck/smMLCK/md/2_mlck_cam
systems=(mcam cam_k30e cam_g40d cam_m36i cam_t34k cam_dN scam1 scam4)

for sys in ${systems[@]}
do

for complex in 2
do

for rep in {1..50}
do

if [ $sys == "scam4" ]
then	
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${complex}_noWAT.top parm [top]
trajin $path/$sys/shortMD/analyses/all_noWAT_${complex}_${rep}.nc 1 -1 parm [top]
rms first :144-404&@CA,C,N,O out
distance dist :269@CG :289@CG out dist_${sys}_${complex}_${rep}.dat
EOD

elif [ $sys == "cam_dN" ]
then
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${complex}_noWAT.top parm [top]
trajin $path/$sys/shortMD/analyses/all_noWAT_${complex}_${rep}.nc 1 -1 parm [top]
rms first :139-399&@CA,C,N,O out
distance dist :264@CG :284@CG out dist_${sys}_${complex}_${rep}.dat
EOD

else
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${complex}_noWAT.top parm [top]
trajin $path/$sys/shortMD/analyses/all_noWAT_${complex}_${rep}.nc 1 -1 parm [top]
rms first :143-403&@CA,C,N,O out
distance dist :268@CG :288@CG out dist_${sys}_${complex}_${rep}.dat
EOD

fi

cpptraj -i tmp.trajin

done
done 
done

rm -rf tmp.trajin
