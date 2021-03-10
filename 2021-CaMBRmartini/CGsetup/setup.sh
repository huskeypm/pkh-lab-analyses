#!/bin/bash

cwd=`pwd`


# generate 200 trials with CaMBR/CaM distance > 40 A
n=1

while [ $n -lt 200 ]
#for i in 1
do
cd $cwd && rm -rf temp  && mkdir temp && cd temp

cp ../input/* .

./construct.sh


#check if the CaMBR is 4.0 nm away from CaM, this is the b-radius at which 
# the interaction between the two binding parters are centrosymmetric. 
# in the BD sims, the two proteins are initially put this far away from each other


# Here we use a more stick criteria, 5.0 nm

rm -rf temp1

echo | awk '{

i=0

camatom = 0
cambratom = 0

while (( getline < "solvtMin.gro") > 0) {
i++
if ( i >= 3 && i <= 302 ) {
        camatom++
        CaM[camatom,0] = $4
        CaM[camatom,1] = $5
        CaM[camatom,2] = $6
 }
if ( i >= 303 && i <= 356 ) {
        cambratom++
        CaMBR[cambratom,0] = $4
        CaMBR[cambratom,1] = $5
        CaMBR[cambratom,2] = $6
 }
 }


Mindis=1000
for(m=1;m<=camatom;m++){
        for(n=1;n<=cambratom;n++) {
        dis = sqrt((CaM[m,0] - CaMBR[n,0])^2 + (CaM[m,1] - CaMBR[n,1])^2 + (CaM[m,2] - CaMBR[n,2])^2)
#       print m,n,dis
        if ( dis <= Mindis ) {
          Mindis= dis
        }
     }
   }


if ( Mindis >= 4.0 ) { 
        print "1" > "temp1"
 } else print "0" > "temp1"

}'

flag=`awk '{ print $1}' temp1`


cat temp1

cd $cwd

if [ $flag -eq 1 ] ; then
	n=$((n+1))
   cp -r temp $n

  cd $cwd/$n
  gmx grompp -f equilibration.mdp -p system.top -c solvtMin.gro -r solvtMin.gro -o equi.tpr -maxwarn 2 2> t7log
  gmx mdrun -deffnm equi -v 2> t8log
else
	rm -r temp
fi
	

done
