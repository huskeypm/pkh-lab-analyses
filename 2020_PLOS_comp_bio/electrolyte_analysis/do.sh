#!/bin/bash

#Get basename from command line
if [ ${#1} = 0 ]
then
  echo "Must provide basename"
  exit 1
fi
basename=$1

nohup /bin/bash ./${basename}.sh >logs/${basename}.txt 2>&1 &