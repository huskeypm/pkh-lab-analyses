#!/usr/bin/python3
"""Read local diffusion coefficient calculation and combine into a single file.

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
import local_D_collection as collector

#Constants
joblist_dtype_cols={'jobname':str,'atm_x':np.float64,'atm_y':np.float64,'atm_z':np.float64}
default_joblist_file="joblist.csv"
default_outdir="collected"
default_dlocal_fname="810_results_D_local.yaml"
default_missing_fname="missing_files.yaml"
default_badfits_fname="bad_fits.yaml"
default_out_fname="Dlocal_collected.csv"

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('setdir',help="Job set directory")
parser.add_argument('--joblist',default=None,help="Path to joblist file (CSV format), defaults to '%s' in setdir"%default_joblist_file)
parser.add_argument('--dlocal_fname',default=default_dlocal_fname,help="Name of the D local results file (yaml format) in each job folder")
parser.add_argument('--outdir',default=None,help="Subdirectory of setdir to use for default output file locations, default is '%s' under setdir"%default_outdir)
parser.add_argument('--missing',nargs="?",default=None,const=False,help="File for storing list of jobs missing the results file, default is '%s' under outdir, provide argument but no file to suppress"%default_missing_fname)
parser.add_argument('--badfits',nargs="?",default=None,const=False,help="File for storing list of bad curve fits, default is '%s' under outdir, provide argument but no file to suppress"%default_badfits_fname)
parser.add_argument('--outfile',default=None,help="Path to the output file (CSV format), default is '%s' under outdir"%default_out_fname)
cmdline=parser.parse_args()

#Check existence of input files and directories
#Job set directory
assert osp.isdir(cmdline.setdir), "Specified job set directory does not exist: %s"%cmdline.setdir
setdir=osp.abspath(cmdline.setdir)
#Job list file
if cmdline.joblist is None:
  jobfile=osp.join(setdir,default_joblist_file)
else:
  jobfile=osp.abspath(cmdline.joblist)
assert osp.isfile(jobfile), "Job list file does not exist: %s"%jobfile
#Output directory
if cmdline.outdir is None:
  outdir=osp.join(setdir,default_outdir)
else:
  outdir=osp.abspath(cmdline.outdir)
outdir_needed=False
for var in [cmdline.missing,cmdline.badfits,cmdline.outfile]:
  outdir_needed = outdir_needed or (var is None)
if outdir_needed:
  if not osp.isdir(outdir):
    os.makedirs(outdir)
  assert osp.isdir(outdir)
#File to list jobs with missing results
if cmdline.missing is False:
  misspath=None
elif cmdline.missing is None:
  misspath=osp.join(outdir,default_missing_fname)
else:
  misspath=osp.abspath(cmdline.missing)
#File to list jobs with failed curve fits
if cmdline.badfits is False:
  badfits=None
elif cmdline.badfits is None:
  badfits=osp.join(outdir,default_badfits_fname)
else:
  badfits=osp.abspath(cmdline.badfits)
#Output file
if cmdline.outfile is None:
  outfpath=osp.join(outdir,default_out_fname)
else:
  outfpath=osp.abspath(cmdline.outfile)

#Load the job list dataframe
jobs_df=pd.read_csv(jobfile,dtype=joblist_dtype_cols)

#Do the collection
output_df,missing_files=collector.do_collection(jobs_df,setdir,cmdline.dlocal_fname)

#If requested, output list of jobs missing the results file
if misspath is not None:
  collector.write_missing(missing_files,misspath)

#If requested, generate list of bad fits
if badfits is not None:
  problem_list=collector.list_bad_fits(output_df)
  collector.write_bad_fits(problem_list,badfits)

#Write the output file
output_df.to_csv(outfpath,index=False)
