#!/usr/bin/python3
"""Compute local diffusion coefficients from trajectories and ACFs.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import sys
import argparse
import os.path as osp
from collections import namedtuple

#Site packages
import numpy as np
import pandas as pd
import matplotlib as mpl
#Set matplotlib backend that won't require Tk
if mpl.get_backend() == 'TkAgg':
  mpl.use("Agg") #this line must appear before the following, and you can't change it later
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import trapz

#Local
mod_dir=osp.abspath(osp.join(osp.split(__file__)[0],"../modules"))
sys.path.append(mod_dir)
import xvg_conversion as xvg

#Constants
basic_columns=['t','x','y','z']
corrheadings=[['t','x'], ['t','y'], ['t','z']]
acf_structure=[2,2,2]
traj_structure=[4]
arg_order=['A','alpha','sigma']
fitparams=namedtuple('fitparams',arg_order)
errortup=fitparams(*[None for itm in arg_order])
errconst="fitting_error"
fitguess=(1,30,1)
rel_diff_tol=0.25
report_template="""%YAML 1.2
---
fit_parameters:
  x:
    A: {fits[x].A}
    alpha: {fits[x].alpha}
    sigma: {fits[x].sigma}
  y:
    A: {fits[y].A}
    alpha: {fits[y].alpha}
    sigma: {fits[z].sigma}
  z:
    A: {fits[z].A}
    alpha: {fits[z].alpha}
    sigma: {fits[z].sigma}
standard_deviations:
  x:
    A: {stdevs[x].A}
    alpha: {stdevs[x].alpha}
    sigma: {stdevs[x].sigma}
  y:
    A: {stdevs[y].A}
    alpha: {stdevs[y].alpha}
    sigma: {stdevs[z].sigma}
  z:
    A: {stdevs[z].A}
    alpha: {stdevs[z].alpha}
    sigma: {stdevs[z].sigma}
msd_values:
  x: {msd[x]}
  y: {msd[y]}
  z: {msd[z]}
fit_results:
  tau:
    x: {tau_fit[x]}
    y: {tau_fit[y]}
    z: {tau_fit[z]}
  D:
    x: {D_fit[x]}
    y: {D_fit[y]}
    z: {D_fit[z]}
numeric_results:
  tau:
    x: {tau_numeric[x]}
    y: {tau_numeric[y]}
    z: {tau_numeric[z]}
  D:
    x: {D_numeric[x]}
    y: {D_numeric[y]}
    z: {D_numeric[z]}
match_ok:
  x: {match_ok[x]}
  y: {match_ok[y]}
  z: {match_ok[z]}
"""

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('acf_xvgfile',help="Path to file containing the acf data (xvg format from gmx analyze")
parser.add_argument('traj_xvgfile',help="Path to file containing the trajectory data (xvg format from gmx traj")
parser.add_argument('outfile',help="Path to the output file (yaml format)")
parser.add_argument('plotfile',help="Path to the output plot file (pdf format)")
cmdline=parser.parse_args()

#Check existence and structure of the xvg files
#(acf)
assert osp.isfile(cmdline.acf_xvgfile), "Specified input file does not exist: %s"%cmdline.acf_xvgfile
found_structure=xvg.get_structure(cmdline.acf_xvgfile)
assert found_structure==acf_structure, "Unexpected xvg file structure: %s"%found_structure
#(trajectory)
assert osp.isfile(cmdline.traj_xvgfile), "Specified input file does not exist: %s"%cmdline.traj_xvgfile
found_structure=xvg.get_structure(cmdline.traj_xvgfile)
assert found_structure==traj_structure, "Unexpected xvg file structure: %s"%found_structure

#Load the data
#(acf)
corrsets=xvg.read_all_sets(cmdline.acf_xvgfile)
corrframes=xvg.sets_to_frames(corrsets,corrheadings)
corrdf=xvg.combine_frames(corrframes,['t'])
#(trajectory)
trajsets=xvg.read_all_sets(cmdline.traj_xvgfile)
trajframes=xvg.sets_to_frames(trajsets,[basic_columns])
intraj=trajframes[0]

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

#Compute MSDs
msdvals={}
for col in 'xyz':
  msdvals[col]=dtraj[col+'sq'].mean()

#Fit to the ACF
def fitfunc(t,A,alpha,sigma):
  return (1.0-A)*np.exp(-alpha*t**2)+A*np.exp(-t/sigma)
results_dict={}
perrs_dict={}
tdata=corrdf['t']
for col in 'xyz':
  try:
    popt,pcov=curve_fit(fitfunc,tdata,corrdf[col],p0=fitguess)
    results_dict[col]=fitparams(*popt)
    perrs_dict[col]=fitparams(*np.sqrt(np.diag(pcov)))
  except RuntimeError:
    results_dict[col]=errortup
    perrs_dict[col]=errortup

#Analytical integrals
tau_analytical={}
Dvals_analytical={}
for p,ptup in results_dict.items():
  if ptup == errortup:
    tau_analytical[p]=errconst
    Dvals_analytical[p]=errconst
  else:
    tau=((1.0-ptup.A)/2.0)*np.sqrt(np.pi/ptup.alpha)+ptup.A*ptup.sigma
    tau_analytical[p]=tau
    Dvals_analytical[p]=msdvals[p]/tau

#Numerical integrals
nstop=1250
tau_numeric={}
Dvals_numeric={}
tstop_dict={}
for col in 'xyz':
  tstop_dict[col]=tdata[nstop]
  tau=trapz(corrdf[col][:nstop],tdata[:nstop])
  tau_numeric[col]=tau
  Dvals_numeric[col]=msdvals[col]/tau

#Look at differences between numerical and analytical results
tau_rel_diffs={}
for col,t_an in tau_analytical.items():
  if t_an == errconst:
    diff=errconst
    tau_rel_diffs[col]=errconst
  else:
    diff=t_an-tau_numeric[col]
    tau_rel_diffs[col]=diff/t_an
rel_diff_ok={}
for col,val in tau_rel_diffs.items():
  if val == errconst:
    rel_diff_ok[col]=False
  else:
    rel_diff_ok[col]=abs(val)<rel_diff_tol

#Output results
reportdict=dict(fits=results_dict,
                stdevs=perrs_dict,
                msd=msdvals,
                tau_fit=tau_analytical,
                D_fit=Dvals_analytical,
                tau_numeric=tau_numeric,
                D_numeric=Dvals_numeric,
                match_ok=rel_diff_ok)
report=report_template.format(**reportdict)
with open(cmdline.outfile,'w') as fp:
  fp.write(report)

#Plots
fig,axlist=plt.subplots(3,1,figsize=(8.5,11))
for ctr,col in enumerate('xyz'):
  ax=axlist[ctr]
  tmax=tstop_dict[col]
  ptup=results_dict[col]
  o=ax.plot(tdata,corrdf[col],'r-',label="data")
  if ptup != errortup:
    fitvals=fitfunc(tdata,*ptup)
    o=ax.plot(tdata,fitvals,'b--',label="fit")
  o=ax.set_xlim(0,tmax)
  o=ax.hlines(0,0,tmax,'k')
  o=ax.legend(loc="upper right")
  o=ax.set_title(col.upper())
  o=ax.set_ylabel("normalized ACF")
  if ctr >= 2:
    o=ax.set_xlabel("time [psec]")
fig.savefig(cmdline.plotfile)