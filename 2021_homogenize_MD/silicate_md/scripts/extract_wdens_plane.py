#!/usr/bin/python3
"""Extract water density data for a single plane aligned with a coordinate axis.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import sys
import argparse
import os
import os.path as osp

#Site packages
import numpy as np
import pandas as pd

#Local
mod_dir=osp.abspath(osp.join(osp.split(__file__)[0],"../modules"))
sys.path.append(mod_dir)
import extract_plane

#Constants

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('infile',help="Path to the input file (OpenDX format)")
parser.add_argument('outfile',help="Path to the output file (pickle format)")
parser.add_argument('coord_index',type=int,help="Axis index (0=x,1=y,2=z) for selected plane")
parser.add_argument('coord_value',type=float,help="Ordinate value for selected plane")
cmdline=parser.parse_args()

#Check that input file exists
assert osp.isfile(cmdline.infile), "Specified input file does not exist: %s"%cmdline.infile

#Perform extraction
extract_plane.top_extraction(cmdline.infile,cmdline.coord_index,cmdline.coord_value,cmdline.outfile)
