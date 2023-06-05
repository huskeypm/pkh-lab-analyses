#!/bin/bash

path=/home/xfang2/backup/faust_backup/mlck/smMLCK/md/2_mlck_cam
systems=(mcam cam_k30e cam_g40d cam_m36i cam_t34k cam_dN scam1 scam4)

for sys in ${systems[@]}
do

for complex in 2
do

for reg in {1..5}
do 

# get hbond
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${complex}_noWAT.top parm [top]
parm $path/$sys/complex2/leap/complex2.top [refTop]
reference /home/xfang2/backup/faust_backup/mlck/smMLCK/md/2_mlck_cam/$sys/complex2/leap/complex2.crd parm [refTop] [crd]
trajin $path/analyses/shortMD/pca/wtAsRef/mlck/cluster_downSampled/${sys}_${complex}_region$reg.nc 1 -1 parm [top]
EOD

if [ $sys == "scam4" ]
then
cat >> tmp.trajin << EOD
nativecontacts name NC1 :1-143&!@H= :144-444&!@H= byresidue out stats_${sys}_${complex}_region$reg.dat mindist maxdist distance 5.0 \\
resout resout_${sys}_${complex}_region$reg.dat writecontacts byAtm_${sys}_${complex}_region$reg.dat \\
reference \\
series seriesout nativeSeries_${sys}_${complex}_region$reg.dat \\
EOD
NUMRES=444

elif [ $sys == "cam_dN" ]
then
cat >> tmp.trajin << EOD
nativecontacts name NC1 :1-138&!@H= :139-439&!@H= byresidue out stats_${sys}_${complex}_region$reg.dat mindist maxdist distance 5.0 \\
resout resout_${sys}_${complex}_region$reg.dat writecontacts byAtm_${sys}_${complex}_region$reg.dat \\
reference \\
series seriesout nativeSeries_${sys}_${complex}_region$reg.dat \\
EOD
NUMRES=439

else
cat >> tmp.trajin << EOD
nativecontacts name NC1 :1-142&!@H= :143-443&!@H= byresidue out stats_${sys}_${complex}_region$reg.dat mindist maxdist distance 5.0 \\
resout resout_${sys}_${complex}_region$reg.dat writecontacts byAtm_${sys}_${complex}_region$reg.dat \\
reference \\
series seriesout nativeSeries_${sys}_${complex}_region$reg.dat \\
EOD
NUMRES=443
fi

cpptraj -i tmp.trajin

sort -n resout_${sys}_${complex}_region$reg.dat | grep -vE "Contacts" > byRes_${sys}_${complex}_region$reg.dat

# post-processing (package data into matrix)
cat > tmp.sh << EOD
echo | awk '{
 while ( ( getline < "byRes_${sys}_${complex}_region$reg.dat" ) > 0 ) {
 printf "%3d %3d %5.2f\n",\$1,\$2,\$3
 printf "%3d %3d %5.2f\n",\$2,\$1,\$3
 }
 }' > dx

 echo | awk '{
  i=0
 while ( ( getline < "dx" ) > 0 ) {
  matrix[1,++i]=\$1
  matrix[2,i]=\$2
  matrix[3,i]=\$3
 }
 N=i

 for(j=1;j<=$NUMRES;j++){    # here \$numRes is the
  for(k=1;k<=$NUMRES;k++){   # total number of residues in the system
       flag=0
       for(z=1;z<=N;z++){
       if( j == matrix[1,z] && k == matrix[2,z] ) {
            flag=1
            v=matrix[3,z]
         }
       }

       if ( flag == 1 ) { print j,k,v
        } else print j,k,"0.00"
 }
  print " "
 }
 }' > pp_${sys}_${complex}_region$reg
EOD

chmod 755 tmp.sh
./tmp.sh

done
done
done

rm -rf tmp*
