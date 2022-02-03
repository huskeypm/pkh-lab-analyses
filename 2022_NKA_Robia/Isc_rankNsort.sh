#!/bin/bash

# print Isc and file name and rank them by Isc
awk '{print $4 " " $NF}' score_dock.sc | sort -n -o IscRanked.dat

# sort by file name
sort -k2 IscRanked.dat | awk '{print $1}' > IscSorted.dat
