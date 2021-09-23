"""For converting wham output into a dataframe with all three coordinates

Expected data structures:

Windows: data about the metafiles that were given to wham.
Each simulation job is 1 window.
For 1D wham windows are grouped into lines, so each metafile contains a line, and
for 2D wham windows are grouped into planes, so each metafile contains a plane.
The "windows" object is a dictionary: {metafile id string: list of jobs}.
Each item in a list of jobs is a sequence of the following: (jobname, coordinate1, ...).
Groups for 1D will have a single coordinate, groups for 2D will have 2 coordinates.

Note that the metafiles prepared for wham are not used here.
Only the metadata (as a yaml file) and wham output files are used.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os.path as osp

#Site packages
import numpy as np
import pandas as pd
from ruamel.yaml import YAML
yaml=YAML(typ="safe",pure=True)
import scipy.constants
physical_constants=scipy.constants.physical_constants

#Constants
outheader=['X','Y','Z','PMF']
coordcols=outheader[0:3]
pmfcol=outheader[3]
selcols_dict={
  1:(0,1,3),
  2:(0,1,2,3)
}
fnametmpl_dict={
  1:"wham_output%s.dat",
  2:"wham2d_output%s.dat"
}
wham_datacols_dict={
  1:["Coor","Free"],
  2:["X","Y","Free"]
}
job_coord_cols=['atm_x','atm_y','atm_z']
dimensionality_mismatch_tmpl="Dimensionality mismatch: %d variable coordinates plus %d constant does not equal %d dimensions"

#Physical constants for unit conversion
#Boltzmann constant
kB_all=physical_constants['Boltzmann constant']
kB=kB_all[0] #J/K
assert kB_all[1]=='J K^-1'
#Avogadro's number
nA_all=physical_constants['Avogadro constant']
nA=nA_all[0]
assert nA_all[1]=='mol^-1'

#For reading the window metadata file
def readfile(fpath,multidoc=False):
  """Read a yaml file into an object"""
  with open(fpath,"r") as fp:
    dat=fp.read()
  if multidoc:
    allobj=yaml.load_all(dat)
  else:
    allobj=yaml.load(dat)
  return allobj

#For reading a single wham output file and converting it to a dataframe (line or plane)
def dat_to_df_arb(fpath,selcols):
  """Read the dat format used by wham"""
  with open(fpath,"r") as fp:
    fdat=fp.readlines()
  datblock=[line.split() for line in fdat if line[0]!="#"]
  header=fdat[0][1:].split()
  columns=[header[sc] for sc in selcols]
  seldata=[[float(row[sc]) for sc in selcols] for row in datblock if len(row)>0]
  df=pd.DataFrame(seldata,columns=columns)
  return df

def dat_to_df_nd(fpath,nd):
  """Call dat_to_df for the specified number of dimensions"""
  return dat_to_df_arb(fpath,selcols_dict[nd])

#For reading all the wham output files for a given set of jobs
def read_all_dats_arb(whamdir,windows,fnametmpl,selcols):
  """Read all the wham output files and return a dictionary of dataframes
  
  Returns a dictionary of dataframes by metafile id."""
  pmf_dfs={}
  for mfid in windows.keys():
    fname=fnametmpl%mfid
    fpath=osp.join(whamdir,fname)
    pmf_dfs[mfid]=dat_to_df_arb(fpath,selcols)
  return pmf_dfs

def read_all_dats_nd(whamdir,windows,nd):
  fnametmpl=fnametmpl_dict[nd]
  selcols=selcols_dict[nd]
  return read_all_dats_arb(whamdir,windows,fnametmpl,selcols)

#For combining wham output into a single dataframe
def combine_pmf_frames(windows,jobs_df,pmf_dfs,nd):
  """Combine pmf dataframes from multiple wham outputs into a single dataframe.

  Arguments:

    - windows = dictionary as described in "Expected data structures" in the module docstring
    - jobs_df = dataframe of jobs, (jobname,atm_x,atm_y,atm_z)
    - pmf_dfs = dictionary of PMF dataframes by metafile id
    - nd = number of dimensions in the wham analysis (1 or 2)
  
  Returns a single dataframe.
  """
  wham_datacols=wham_datacols_dict[nd]
  rows=[]
  #Loop over metafiles
  for mfid,joblist in windows.items():
    #Get the subset of the jobs that were part of this metafile
    jobnames=[job[0] for job in joblist]
    jobsubdf=jobs_df.loc[jobs_df['jobname'].isin(jobnames)]
    #Get the coordinates that are constant for this metafile
    constpart=[]
    for col in job_coord_cols:
      uniquevals=jobsubdf[col].unique()
      if uniquevals.shape[0]==1:
        #This coordinate is constant
        constpart.append(uniquevals[0])
    assert nd+len(constpart)==len(job_coord_cols), dimensionality_mismatch_tmpl%(nd,len(constpart),len(job_coord_cols))
    #Add all the points from this metafile
    for datatup in pmf_dfs[mfid].itertuples():
      varpart=[getattr(datatup,col) for col in wham_datacols]
      rows.append(constpart+varpart)
  pmf=pd.DataFrame(rows,columns=outheader)
  return pmf

#For unit conversion
def convert_to_kBT(df,temp):
  """Convert from units of kJ/mol to kBT

  Arguments:

    - df = collected dataframe of PMF data
    - temp = temperature, as number, in Kelvin

  Returns:

    - out_df = dataframe with the PMF column converted to units of kBT."""
  kBT = kB * temp / 1000 # kJ
  to_kBT = 1.0/kBT * 1.0/nA #kBT / (kJ/mol)
  outdata=[]
  for row in df.itertuples(index=False,name="pmfdata"):
    outrow=[getattr(row,k) for k in coordcols]
    outrow.append(getattr(row,pmfcol)*to_kBT)
    outdata.append(outrow)
  out_df=pd.DataFrame(outdata,columns=outheader)
  return out_df
