#!/bin/bash

path=/home/xfang2/backup/faust_backup/mlck/smMLCK/md/2_mlck_cam
systems=(scam4 mcam cam_k30e cam_g40d cam_dN cam_m36i cam_t34k scam1)
nRep=50
start=5
end=443
covarFile=WT_covar.dat

for sys in ${systems[@]}
do

for comp in 2
do	

# 1. get WT covar matrix (only CaM)
if [ $sys == mcam ]
then
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${comp}_noWAT.top parm [top]
parm $path/$sys/complex2/leap/complex2.top [refTop]
reference $path/$sys/complex2/leap/complex2.crd parm [refTop] [crd]
trajin $path/$sys/shortMD/analyses/all_${comp}.nc 1 -1 parm [top]
EOD

cat >> tmp.trajin << EOD
rms test ref [crd] :143-403@CA,C,N,O
matrix covar name covar :${start}-${end}@CA,C,N,O
diagmatrix covar out $covarFile vecs 20 name myEvecs nmwiz nmwizvecs 3 nmwizfile dna.nmd nmwizmask :${start}-${end}@CA,C,N,O
EOD

cpptraj -i tmp.trajin

# reformat PCA weights file
rm -rf temp
rm -rf PCA_weight PCA_weight_tmp
grep -n '\*' $covarFile | awk -F ":" '{ print $1 + 1 }' >> temp

n=`cat temp | wc -l`

for i in `seq 1 $n`
do
line_number=`awk '{ if ( NR == '$i' ) print $1 }' temp`
awk '{ if ( NR == '$line_number' ) print '$i',$2 }' $covarFile >> PCA_weight_tmp
done

sum=`awk 'BGEIN{ sum = 0 } { sum += $2 }END{print sum}' PCA_weight_tmp`

awk '{ print $1,$2/'$sum'}' PCA_weight_tmp > PCA_weight

rm -rf tmp.* *temp*

fi

# 2. project traj onto cov matrix
cat > tmp.trajin << EOD
parm $path/$sys/shortMD/analyses/${comp}_noWAT.top parm [top]
parm $path/$sys/complex2/leap/complex2.top [refTop]
reference $path/$sys/complex2/leap/complex2.crd parm [refTop] [crd]
trajin $path/$sys/shortMD/analyses/all_${comp}.nc 1 -1 parm [top]
EOD

if [ $sys == "scam4" ]
then
cat >> tmp.trajin << EOD
readdata $covarFile
rms test ref [crd] :144-404@CA,C,N,O
projection md_project evecs $covarFile out md_project_${sys}_${comp} beg 1 end 3 :$((start+1))-$((end+1))@CA,C,N,O

#1D histograms (open with xmgrace)
hist md_project:1 free 310 bins 50 out ${sys}_${comp}_hists1.agr name PC1
hist md_project:2 free 310 bins 50 out ${sys}_${comp}_hists2.agr name PC2
hist md_project:3 free 310 bins 50 out ${sys}_${comp}_hists3.agr name PC3

#2D histograms (open with gnuplot)
hist md_project:1 md_project:2 bins 50 out ${sys}_${comp}_hists_1-2.gnu name PC12 free 310
hist md_project:1 md_project:3 bins 50 out ${sys}_${comp}_hists_1-3.gnu name PC13 free 310
hist md_project:2 md_project:3 bins 50 out ${sys}_${comp}_hists_2-3.gnu name PC23 free 310
EOD

elif [ $sys == "cam_dN" ]
then
cat >> tmp.trajin << EOD
readdata $covarFile
rms test ref [crd] :139-399@CA,C,N,O
projection md_project evecs $covarFile out md_project_${sys}_${comp} beg 1 end 3 :1-$((end - 4))@CA,C,N,O

#1D histograms (open with xmgrace)
hist md_project:1 free 310 bins 50 out ${sys}_${comp}_hists1.agr name PC1
hist md_project:2 free 310 bins 50 out ${sys}_${comp}_hists2.agr name PC2
hist md_project:3 free 310 bins 50 out ${sys}_${comp}_hists3.agr name PC3

#2D histograms (open with gnuplot)
hist md_project:1 md_project:2 bins 50 out ${sys}_${comp}_hists_1-2.gnu name PC12 free 310
hist md_project:1 md_project:3 bins 50 out ${sys}_${comp}_hists_1-3.gnu name PC13 free 310
hist md_project:2 md_project:3 bins 50 out ${sys}_${comp}_hists_2-3.gnu name PC23 free 310
EOD

else
cat >> tmp.trajin << EOD
readdata $covarFile
rms test ref [crd] :143-403@CA,C,N,O
projection md_project evecs $covarFile out md_project_${sys}_${comp} beg 1 end 3 :${start}-${end}@CA,C,N,O

#1D histograms (open with xmgrace)
hist md_project:1 free 310 bins 50 out ${sys}_${comp}_hists1.agr name PC1
hist md_project:2 free 310 bins 50 out ${sys}_${comp}_hists2.agr name PC2
hist md_project:3 free 310 bins 50 out ${sys}_${comp}_hists3.agr name PC3

#2D histograms (open with gnuplot)
hist md_project:1 md_project:2 bins 50 out ${sys}_${comp}_hists_1-2.gnu name PC12 free 310
hist md_project:1 md_project:3 bins 50 out ${sys}_${comp}_hists_1-3.gnu name PC13 free 310
hist md_project:2 md_project:3 bins 50 out ${sys}_${comp}_hists_2-3.gnu name PC23 free 310
hist md_project:1 md_project:2 bins 50 out ${sys}_${comp}_hists_prob_1-2.gnu name PC12p norm
hist md_project:1 md_project:3 bins 50 out ${sys}_${comp}_hists_prob_1-3.gnu name PC13p norm
hist md_project:2 md_project:3 bins 50 out ${sys}_${comp}_hists_prob_2-3.gnu name PC23p norm
EOD
fi

cpptraj -i tmp.trajin

done
done
