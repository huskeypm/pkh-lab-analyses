#!/usr/bin/python3
"""Compute deviations from a trajectory xvg file, and optionally compute squares and ACFs.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import sys
import argparse
import os.path as osp

#Site packages
import numpy as np
import pandas as pd

#Local
mod_dir=osp.abspath(osp.join(osp.split(__file__)[0],"../modules"))
sys.path.append(mod_dir)
import xvg_conversion as xvg

#Constants
basic_columns=['t','x','y','z']

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('xvgfile',help="Path to file containing the trajectory data (xvg format from gmx traj")
parser.add_argument('outfile',help="Path to the output file (CSV format)")
parser.add_argument('--shift',choices=['mean','initial','none'],default='mean',help="Take deviation from mean position (default), initial position, or use raw coordinate")
parser.add_argument('--calc_square',action="store_true",help="To also calculate squared deviations")
parser.add_argument('--calc_acf',action="store_true",help="To also perform calculation of ACF (note: this can be time-consuming, and the results at large times may not agree with other methods)")
cmdline=parser.parse_args()

#Check existence and structure of the xvg file
assert osp.isfile(cmdline.xvgfile), "Specified input file does not exist: %s"%cmdline.xvgfile
found_structure=xvg.get_structure(cmdline.xvgfile)
assert found_structure==[4], "Unexpected xvg file structure: %s"%found_structure

#Load the data
trajsets=xvg.read_all_sets(cmdline.xvgfile)
trajframes=xvg.sets_to_frames(trajsets,[basic_columns])
intraj=trajframes[0]

#Compute the deltas from the mean position
dtraj=pd.DataFrame(columns=basic_columns)
dtraj['t']=intraj['t']
for col in 'xyz':
  if cmdline.shift == 'mean':
    shift=intraj[col].mean()
  elif cmdline.shift == 'initial':
    shift=intraj[col][0]
  elif cmdline.shift == 'none':
    shift=0
  dtraj[col]=intraj[col]-shift

#Compute the squared deviation versus time
def dosquare(row,col):
  return row[col]**2

if cmdline.calc_square:
  for col in 'xyz':
    sqf = lambda row: dosquare(row,col)
    dtraj[col+'sq']=dtraj.apply(sqf,axis=1) #axis=0 to pass columns, axis=1 to pass rows

#ACF calculation
if cmdline.calc_acf:
  for c in 'xyz':
    col='ACF'+c
    dtraj[col]=np.nan
    for idx in range(len(dtraj)):
    #for idx in range(4000):
      dtraj.loc[idx,col]=dtraj[c].autocorr(idx)

#Write the output file
dtraj.to_csv(cmdline.outfile,index=False)
