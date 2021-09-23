#!/usr/bin/python3
"""Compute square deviations from a trajectory xvg file

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
corrheadings=[['t','x'], ['t','y'], ['t','z']]

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('xvgfile',help="Path to file containing the acf data (xvg format from gmx analyze")
parser.add_argument('outfile',help="Path to the output file (CSV format)")
cmdline=parser.parse_args()

#Check existence and structure of the xvg file
assert osp.isfile(cmdline.xvgfile), "Specified input file does not exist: %s"%cmdline.xvgfile
found_structure=xvg.get_structure(cmdline.xvgfile)
assert found_structure==[2,2,2], "Unexpected xvg file structure: %s"%found_structure

#Load the data
corrsets=xvg.read_all_sets(cmdline.xvgfile)
corrframes=xvg.sets_to_frames(corrsets,corrheadings)
corrdf=xvg.combine_frames(corrframes,['t'])

#Write the output file
corrdf.to_csv(cmdline.outfile,index=False)
