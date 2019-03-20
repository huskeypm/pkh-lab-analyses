#!/bin/bash

 input=$1
 echo | awk '{

 while ( ( getline < "'$input'" ) > 0 ) {
 avge = ($2+$3 + $4)/3
 std = sqrt(((avge - $2)^2 + (avge -$3)^2 + (avge-$4)^2)/3)
 print $1,avge, std
 }
 }' 
