#!/bin/bash

path=/home/xfang2/backup/faust_backup/mlck/smMLCK/md/2_mlck_cam
systems=(cam_t34k cam_dN mcam cam_k30e cam_g40d cam_m36i scam1 scam4)

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
distance dist :1-143@CA,C,N,O :144-404@CA,C,N,O out dist_${sys}_${complex}_${rep}.dat
EOD

elif [ $sys == "cam_dN" ]
then
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${complex}_noWAT.top parm [top]
trajin $path/$sys/shortMD/analyses/all_noWAT_${complex}_${rep}.nc 1 -1 parm [top]
rms first :139-399&@CA,C,N,O out
distance dist :1-138@CA,C,N,O :139-399@CA,C,N,O out dist_${sys}_${complex}_${rep}.dat
EOD

else
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${complex}_noWAT.top parm [top]
trajin $path/$sys/shortMD/analyses/all_noWAT_${complex}_${rep}.nc 1 -1 parm [top]
rms first :143-403&@CA,C,N,O out
distance dist :1-142@CA,C,N,O :143-403@CA,C,N,O out dist_${sys}_${complex}_${rep}.dat
EOD

fi

cpptraj -i tmp.trajin

init=$(awk 'NR==2 {print $2}' dist_${sys}_${complex}_${rep}.dat)
awk -v var="$init" 'NR>1 {$(NF+1)=$2 - var;}1' OFS=" " dist_${sys}_${complex}_${rep}.dat > tmp
mv tmp dist_${sys}_${complex}_${rep}.dat


done
done 
done

rm -rf tmp.trajin
