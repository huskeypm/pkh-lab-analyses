"""For collecting local diffusion coefficient data from multiple MD jobs into one file.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os.path as osp
from collections import OrderedDict as odict

#Site packages
import numpy as np
import pandas as pd

#Local packages
import simproc
import simproc.requesthandler.yaml_manager as yaml_manager
import simproc.requesthandler.nested as nested

#Constants
job_data_cols=["job","X","Y","Z"]
dlocal_fname_default="810_results_D_local.yaml"

#Construct column to key sequence mapping
column_paths=odict()
for coord in "xyz":
  column_paths["msd_%s"%coord]=["msd_values",coord]
for rtype in ["fit","numeric"]:
  for vtype in ["D","tau"]:
    for coord in "xyz":
      column_paths["%s_%s_%s"%(vtype,rtype,coord)]=["%s_results"%rtype,vtype,coord]
for coord in "xyz":
  column_paths["ok_%s"%coord]=["match_ok",coord]
start_columns=job_data_cols+list(column_paths.keys())
all_columns=start_columns+["all_ok"]

#Functions for processing a single table row
def get_row(infpath):
  resdict=yaml_manager.readfile(infpath)
  return [nested.get_nested(resdict,dpath) for dpath in column_paths.values()]

def do_check_all(row):
  res=True
  for coord in "xyz":
    res=res and row["ok_%s"%coord]
  return res

#Function for doing the whole collection
def do_collection(jobsdf,setdir,dlocal_fname=dlocal_fname_default):
  """Read all the job result files and combine

  Arguments:

    - jobsdf = job list dataframe
    - setdir = path to directory for the job set
    - dlocal_fname = name of the results yaml file within each job folder in the setdir
  
  Return values:

    - data_df = dataframe of the results
    - missing_files = list of jobs where the results file is missing"""
  #Initializations
  missing_files=[]
  datarows=[]
  #For each job in the set, read the data and put it in the proper columns
  for jobrow in jobsdf.itertuples(index=False,name="jobdata"):
    jobname=jobrow.jobname
    jobcoords=[getattr(jobrow,"atm_%s"%coord) for coord in "xyz"]
    dlocal_fpath=osp.join(setdir,jobname,dlocal_fname)
    if osp.isfile(dlocal_fpath):
      resdata=get_row(dlocal_fpath)
      row=[jobname]+jobcoords+resdata
      datarows.append(row)
    else:
      missing_files.append(jobname)
  #Create dataframe
  data_df=pd.DataFrame(datarows,columns=start_columns)
  #Add column for overall fitting result
  data_df['all_ok']=data_df.apply(do_check_all,axis=1)
  #Done
  return data_df,missing_files

#Provide a list of all the fitting problems
def list_bad_fits(data_df):
  """List all the fitting problems

  Arguments:

    - data_df = dataframe as returned by do_collection
  
  Return value:

    - problem_list = """
  #Select rows where at least one fitting failed
  badrows=data_df.loc[data_df["all_ok"]==False,:]
  #For each row with a failed fitting, list which ones
  problem_list=odict()
  for row in badrows.itertuples(index=False,name="datarow"):
    failures=[coord for coord in "xyz" if not getattr(row,"ok_%s"%coord)]
    problem_list[row.job]=failures
  #Done
  return problem_list

#For writing out calculations to a properly formatted file that is still parseable yaml

def write_bad_fits(problem_list,fpath):
  """Write the problem list in a nice yaml format"""
  with open(fpath,"w") as fp:
    for jobname,failures in problem_list.items():
      outlist_data=",".join(failures)
      outline="%s: [%s]\n"%(jobname,outlist_data)
      fp.write(outline)
  return

def write_missing(missing_files,fpath):
  """Write the list of missing files to a file"""
  with open(fpath,"w") as fp:
    for jobname in missing_files:
      fp.write("- %s\n"%jobname)
  return