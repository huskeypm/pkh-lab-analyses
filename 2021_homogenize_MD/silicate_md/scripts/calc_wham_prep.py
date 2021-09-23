#!/usr/bin/python3
"""Prepare simulation results for WHAM calculations.

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
traj_structure=[4]
headerlines={False: "#t (ps), z (nm)\n",
             True: "#t (ps), y (nm), z (nm)\n"}

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('traj_xvgfile',help="Path to file containing the trajectory data (xvg format from gmx traj")
parser.add_argument('outfile',help="Path to the output file (dat format)")
parser.add_argument('--2D',dest="is2D",action="store_true",help="To generate output for wham-2d instead")
cmdline=parser.parse_args()

#Check existence and structure of the xvg file
#(trajectory)
assert osp.isfile(cmdline.traj_xvgfile), "Specified input file does not exist: %s"%cmdline.traj_xvgfile
found_structure=xvg.get_structure(cmdline.traj_xvgfile)
assert found_structure==traj_structure, "Unexpected xvg file structure: %s"%found_structure

#Load the data
#(trajectory)
trajsets=xvg.read_all_sets(cmdline.traj_xvgfile)
trajframes=xvg.sets_to_frames(trajsets,[basic_columns])
intraj=trajframes[0]

#Create output file
with open(cmdline.outfile,'w') as fp:
  fp.write(headerlines[cmdline.is2D])
  for row in intraj.itertuples():
    if cmdline.is2D:
      fp.write("{:10.3f} {:10.5f} {:10.5f}\n".format(row.t,row.y,row.z))
    else:
      fp.write("{:10.3f} {:10.5f}\n".format(row.t,row.z))
