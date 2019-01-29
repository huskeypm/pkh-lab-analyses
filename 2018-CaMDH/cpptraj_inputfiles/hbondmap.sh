#!/bin/bash

# Post-processing the cpptraj hbond output results
# the input file is "avout" which is specified in the hbond.trajin
# the output file is "pp", which can be loaded in the hbond.ipynb 


# there is an example "avout" file in this directory

# Usage: ./hbondmap.sh



cwd=`pwd`


#----------
# remove the hhbonds fromed between neighbouring residues,
# namely, we only consider hbonds formed betwee i and i+2, i+3.. residues
#------
 awk '{ if ( NR > 1  ) print }' avout > t

 sed -i 's/@/_/g' t 

 awk -F "_" '{ printf "%3d %3d %6.2f\n",$2,$4,substr($7,18,10) }' t > tx

 echo | awk 'function abs(v) {return v < 0 ? -v : v}
{
  while ( ( getline < "tx" ) > 0 ) {
  a=$1-$2
   if ( abs(a) >= 2 && $3 != 0) print  # residue seperation must be > 2
 }
}' > txx




#--------------------------------------
# there may exist multiple hbonds between one residue pair
# , we only keep the hbond with highest life time percentage
# and remove other hbonds, thus making sure only one
# hbond is in each residue pair
#-------------------------------------
 sort -k 1 -g txx > txxx
 
 echo | awk '{
 i=0
 while ( ( getline < "txxx" ) > 0 ) {
  fid[++i]=$1
  sid[i]=$2
  v[i]=$3
 if ( i > 1 ) {
 if ( fid[i] == fid[i-1] && sid[i] == sid[i-1] ) {
 flag=0
 } else flag=1

 if ( flag ==1 ) print fid[i-1],sid[i-1],v[i-1]
 }
 }
 }' > ax
 
   echo | awk '{
  i=0
 while ( ( getline < "ax" ) > 0 ) {
  fid[++i]=$1
  sid[i]=$2
  v[i]=$3
 }
 N=i

 for(j=1;j<=N;j++){
  flag=0
  for(k=1;k<=N;k++){
       if( fid[j] == sid[k] && sid[j] == fid[k] ) {
            flag=1
          if ( v[j] >= v[k] ) {
            aa = v[j]
          } else aa = v[k]
          }
       }
     if ( flag == 1 ) {
         if ( fid[j] <= sid[j] ) {
         print fid[j],sid[j],aa
           } else print sid[j],fid[j],aa
       } else print fid[j],sid[j],v[j]
  }
 }' > bxx

 sort -k 1 -g bxx | uniq > cx

#--------------------------
# Lastly, generate a matrix data
# that can be directly loaded by 
# np.loadtext in ipython notebook
#--------------------------


  echo | awk '{
 while ( ( getline < "cx" ) > 0 ) {
 printf "%3d %3d %5.2f\n",$1,$2,$3
 printf "%3d %3d %5.2f\n",$2,$1,$3
 }
 }' > dx
 
 echo | awk '{
  i=0
 while ( ( getline < "dx" ) > 0 ) {
  matrix[1,++i]=$1
  matrix[2,i]=$2
  matrix[3,i]=$3
 }
 N=i
  
 for(j=1;j<=213;j++){    # here "213" is the 
  for(k=1;k<=213;k++){   # total number of residues in the system
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
 }' > pp
