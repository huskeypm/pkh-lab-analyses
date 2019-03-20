#!/bin/bash

rm -rf temp
input1=$1
input2=$2
input3=$3

touch temp

cat $input1 >> temp
cat $input2 >> temp
cat $input3 >> temp

sed -i '/#/d' temp

echo | awk '{ 

i=0
sum=0
while ( ( getline < "temp" ) > 0 ) {
 data[++i]=$2
 sum += $2
 }
 avg = sum/i

 std = 0
 for(k=1;k<=i;k++) {
  std += (data[k] - avg)^2
 }

 print avg,sqrt(std/i)
 }'
