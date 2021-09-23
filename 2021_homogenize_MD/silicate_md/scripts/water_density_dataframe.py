#!/usr/bin/python3
"""Convert water density dx file to a pandas dataframe using grid cell center-points.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import sys
import argparse
import os.path as osp

#Site packages
import numpy as np
import pandas as pd
from gridData import Grid

#Local

#Constants
header=["x","y","z","density"]

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('dxfile',help="Path to file containing the water density data (dx format from gmx traj")
parser.add_argument('outfile',help="Path to the output file (CSV format)")
cmdline=parser.parse_args()

#Check existence of the input file
assert osp.isfile(cmdline.dxfile), "Specified input file does not exist: %s"%cmdline.dxfile

#Load the data
dn=Grid(cmdline.dxfile)

#Convert to dataframe using cell midpoints
def makeframe(dn):
  data=[]
  for idx in np.ndindex(dn.grid.shape):
    cen=dn.delta*np.array(idx)+dn.origin
    val=dn.grid[idx]
    row=np.append(cen,[val])
    data.append(row)
  df=pd.DataFrame(data,columns=header)
  return df

df=makeframe(dn)

#Output
df.to_csv(cmdline.outfile,index=False)
