#!/bin/bash

# Get the contract trace between CaMBR and the N/C domains of CaM

# There will be 3 output files: a) "contact" --records the total number of contacts 
#                               b) "Ncontact"/"Ccontact" -- records the residue pairs forming the contact




cwd=`pwd`


for z in 0M 0.15M 0.5M
do
for i in `seq 101 200`
do
cd $cwd/$z/$i

if [ -f contacts ];then
	continue
else

rm -rf contacts Ncontact Ccontact

echo | awk '{
dt = 1.5 # ns
cutoff = 5.5 # distance cutoff of contact 


Nframe = 0
while ( ( getline < "viz.pdb" ) > 0 ) {
 if ( $1 == "MODEL" ) {
	Nframe++
	ni = 0
	nc = 0
	ncambr = 0
	}
	
 if ( $1 == "ATOM" && $5 == "A" ) {
	if ( $6 <= 71 ) {
        ni++
	N[Nframe,ni,0] = $7
	N[Nframe,ni,1] = $8
	N[Nframe,ni,2] = $9
	N[Nframe,ni,4] = $6 # the residue Number
	}
	if ( $6 > 71 && $6 <= 143 ) {
	nc++
	C[Nframe,nc,0] = $7
	C[Nframe,nc,1] = $8
	C[Nframe,nc,2] = $9
	C[Nframe,nc,4] = $6
	}
	}

 if ( $1 == "ATOM" && $5 == "B" ) {
	ncambr++
	CamBR[Nframe, ncambr,0] = $7
	CamBR[Nframe, ncambr,1] = $8
	CamBR[Nframe, ncambr,2] = $9
	CamBR[Nframe, ncambr,4] = $6 # the residue number
}
}


totalFrame=Nframe
Natom = ni
Catom = nc
cambratom = ncambr

for(m=1;m<=Nframe;m++){
	time = dt*(m-1)
	Ncontacts = 0
	Ccontacts = 0
	# N domain contacts 
	pairN = ""
	pairC = ""
	for(j=1;j<=Natom;j++){
		for(z=1;z<=cambratom;z++){
		dist = sqrt((N[m,j,0] - CamBR[m,z,0])^2 + (N[m,j,1] - CamBR[m,z,1])^2 + (N[m,j,2] - CamBR[m,z,2])^2 )
		if (dist <= cutoff ) {Ncontacts++; pairN = pairN N[m,j,4] "-" CamBR[m,z,4] "  " }
		}
 	}

	for(k=1;k<=Catom;k++){
		for(z=1;z<=cambratom;z++){
		dist = sqrt((C[m,k,0] - CamBR[m,z,0])^2 + (C[m,k,1] - CamBR[m,z,1])^2 + (C[m,k,2] - CamBR[m,z,2])^2 )
		if (dist <= cutoff ) {Ccontacts++;  pairC = pairC C[m,k,4] "-" CamBR[m,z,4] "  "}
		}
	}
        printf "%8.2f %4d %4d\n",time,Ncontacts,Ccontacts
        printf "%8.2f %s\n",time,pairN >> "Ncontact"
        printf "%8.2f %s\n",time,pairC >> "Ccontact"
 }
 }' > contacts
 fi
done
done
