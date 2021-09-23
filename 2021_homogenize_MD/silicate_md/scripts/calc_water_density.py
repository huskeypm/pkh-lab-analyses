#!/usr/bin/python3
"""Compute water density and save to file

The topology and trajectory files must be in formats readable by MDAnalysis.

The list of topology formats is at
https://www.mdanalysis.org/docs/documentation_pages/topology/init.html#supported-topology-formats

The list of trajectory formats is at
https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#id2

A list of valid strings for density units can be found at
https://www.mdanalysis.org/docs/documentation_pages/analysis/density.html#MDAnalysis.analysis.density.Density.convert_density

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import sys
import argparse
import os.path as osp

#Site packages
import numpy as np
import MDAnalysis
from MDAnalysis.analysis.density import DensityAnalysis
from MDAnalysis import units

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument('topology',help="Path to the (compiled) gromacs topology (tpr) file")
parser.add_argument('trajectory',help="Path to the trajectory file (see help for acceptable formats")
parser.add_argument('outfile',help="Path to the output file (OpenDX format)")
parser.add_argument('delta',help="Cell size (size of histogram bins)",nargs="?",default="1.0")
parser.add_argument('units',help="Density units to be used (see help for acceptable values)",nargs="?",default="SPC")
parser.add_argument('--multitraj',nargs="+",help="Additional trajectories to load")
cmdline=parser.parse_args()
delta=float(cmdline.delta)

#If more than one trajectory, generate list
trajlist=[cmdline.trajectory]
if cmdline.multitraj is not None:
  trajlist+=cmdline.multitraj

#Load data
universe = MDAnalysis.Universe(cmdline.topology,trajlist)

#Calculate density
ow = universe.select_atoms("name OW")
da = DensityAnalysis(ow, delta=delta)
da.run()
dn = da.density
dn.convert_density(cmdline.units)

#Export result
dn.export(cmdline.outfile,"dx")
