#!/bin/bash

aa=$1
bb=$2
cc=$3

paste $aa $bb $cc | awk '{ if (NR > 1) {

avge = ($2+$4+$6)/3
std = sqrt((avge - $2)^2 + (avge-$4)^2 + (avge - $6)^2/3)

print $1,avge,std
}
}'
