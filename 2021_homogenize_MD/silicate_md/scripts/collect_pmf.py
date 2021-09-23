#!/usr/bin/python3
"""Read wham output and compile into dataframe, and convert units.

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
import wham_data

#Constants
joblist_dtype_cols={'jobname':str,'atm_x':np.float64,'atm_y':np.float64,'atm_z':np.float64}
default_whamdir_rel="wham"
default_joblist_file="joblist.csv"
windfile_defaults={1:"windows.yaml",2:"windows2d.yaml"}

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('setdir',help="Job set directory")
parser.add_argument('num_dim',type=int,help="Number of dimensions, either 1 or 2")
parser.add_argument('--whamdir',default=None,help="Path to wham data directory, defaults to '%s' in setdir"%default_whamdir_rel)
parser.add_argument('--joblist',default=None,help="Path to joblist file (CSV format), defaults to '%s' in setdir"%default_joblist_file)
parser.add_argument('--windowsfile',default=None,help="Path to windows metadata yaml file, default varies with number of dimensions")
parser.add_argument('--temperature',default=298.0,help="Temperature to use for unit conversion to kBT, in Kelvin, as float, default=298")
parser.add_argument('--outfile',default=None,help="Path to the output file (CSV format), default is in whamdir, with name depending on number of dimensions")
parser.add_argument('--no_conversion',action="store_false",dest="convert",help="To skip the unit conversion step")
cmdline=parser.parse_args()

#Check existence of input files and directories
#Job set directory
assert osp.isdir(cmdline.setdir), "Specified job set directory does not exist: %s"%cmdline.setdir
setdir=osp.abspath(cmdline.setdir)
#Wham directory
if cmdline.whamdir is None:
  whamdir=osp.join(setdir,default_whamdir_rel)
else:
  whamdir=osp.abspath(cmdline.whamdir)
assert osp.isdir(whamdir), "Directory for wham files does not exist: %s"%whamdir
#Job list file
if cmdline.joblist is None:
  jobfile=osp.join(setdir,default_joblist_file)
else:
  jobfile=osp.abspath(cmdline.joblist)
assert osp.isfile(jobfile), "Job list file does not exist: %s"%jobfile
#windows metdata file
if cmdline.windowsfile is None:
  windfile=osp.join(whamdir,windfile_defaults[cmdline.num_dim])
else:
  windfile=osp.abspath(cmdline.windowsfile)
assert osp.isfile(windfile), "Windows yaml file does not exist: %s"%windfile
#output file
if cmdline.outfile is None:
  outfpath=osp.join(whamdir,"pmf%dd.csv"%(cmdline.num_dim))
else:
  outfpath=osp.abspath(cmdline.outfile)

#Load the data
#Job list dataframe
jobs_df=pd.read_csv(jobfile,dtype=joblist_dtype_cols)
#Windows metadata
windows=wham_data.readfile(windfile)
#All wham output data files
pmf_dfs=wham_data.read_all_dats_nd(whamdir,windows,cmdline.num_dim)

#Combine wham output into single dataframe
pmf=wham_data.combine_pmf_frames(windows,jobs_df,pmf_dfs,cmdline.num_dim)

#Do the unit conversion
if cmdline.convert:
  pmf=wham_data.convert_to_kBT(pmf,cmdline.temperature)

#Write the output file
pmf.to_csv(outfpath,index=False)
