#!/bin/bash

grep -v "#" lowd.gmds | awk '{ print $1,$2 }' > temp_project

dimproj -D 32 -d 2 -w -grid 400,20,200 -fun-hd 18,8,6 -fun-ld 18,1,2 -gopt 5 -P landmark -p temp_project < hddata > FULL_PROJECTION 
