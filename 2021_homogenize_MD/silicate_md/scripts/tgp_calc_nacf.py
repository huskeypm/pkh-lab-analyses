#!/usr/bin/python3
"""Compute normalized positional ACFs from a trajectory .dat file

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse

#Site packages
import numpy as np
import pandas as pd

#Constants
basic_columns=['t','x','y','z']

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('datfile',help="Path to file containing the trajectory data: space-delimited t,x,y,z")
parser.add_argument('outfile',help="Path to the output file (CSV format)")
cmdline=parser.parse_args()

#Load the data
intraj=pd.read_csv(cmdline.datfile,sep=" ",names=basic_columns)

#Compute the deltas from the mean position
dtraj=pd.DataFrame(columns=basic_columns)
dtraj['t']=intraj['t']
for col in 'xyz':
  dtraj[col]=intraj[col]-intraj[col].mean()

#Compute the squared deviation versus time
def dosquare(row,col):
  return row[col]**2

for col in 'xyz':
  sqf = lambda row: dosquare(row,col)
  dtraj[col+'sq']=dtraj.apply(sqf,axis=1) #axis=0 to pass columns, axis=1 to pass rows

#ACF calculation
for c in 'xyz':
  col='ACF'+c
  dtraj[col]=np.nan
  for idx in range(len(dtraj)):
  #for idx in range(4000):
    dtraj.loc[idx,col]=dtraj[c].autocorr(idx)

#Write the output file
dtraj.to_csv(cmdline.outfile,index=False)
