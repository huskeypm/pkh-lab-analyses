"""For extracting a single plane of water density data from a dx file.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages
import numpy as np

#Local packages
from gridData import Grid
import simproc
import simproc.requesthandler.pickle_manager as pickle_manager

#Constants
match_errstr_tmpl="Failed coordinate match: %0.3f<%0.3f<%0.3f is not True"

#Functions
def load_dx_convert(infpath):
  """Load the data and convert edges from Angstroms to nm"""
  dn=Grid(infpath)
  grid=dn.grid
  raw_edges=dn.edges
  edges=[ed/10.0 for ed in raw_edges]
  return grid,edges

def extract_plane(in_grid,in_edges,plane_coord_idx,plane_coord_val):
  """Extract a single plane of data from the 3D data.

  Arguments:

    - in_grid = 3D numpy array of data
    - in_edges = sequence of 3 1D arrays for the cell edges
    - plane_coord_idx: index of the cartesian axis perpendicular to the plane
        (that is, the cartesian value that is constant in the plane).
        0 for X-plane, 1 for Y-plane, 2 for Z-plane.
    - plane_coord_val: value for the ordinate of the plane, as float
  
  Returns:

    - out_grid = 2D numpy array of data
    - out_edges = sequence of 2 1D arrays for the cell edges"""
  #Get the index of the plane
  plane_bounds=in_edges[plane_coord_idx]
  bools=[v<plane_coord_val for v in plane_bounds]
  match_index=bools.index(False)
  #Confirm that the search for the plane's index worked
  lower_bound=plane_bounds[match_index-1]
  upper_bound=plane_bounds[match_index]
  match_ok=lower_bound<plane_coord_val and upper_bound>plane_coord_val
  match_errstr=match_errstr_tmpl%(lower_bound,plane_coord_val,upper_bound)
  assert match_ok,match_errstr
  #Extract the data
  slices=[slice(None)]*3
  slices[plane_coord_idx]=match_index
  slices=tuple(slices)
  out_grid=in_grid[slices]
  out_edges=[in_edges[idx] for idx in range(3) if idx!=plane_coord_idx]
  #Done
  return out_grid,out_edges

def save_pickle(grid,edges,outfpath):
  """Save the results to a pickle file"""
  out_dict={
  'grid': grid,
  'edges': edges
  }
  pickle_manager.writefile(out_dict,outfpath)
  return

def top_extraction(infpath,plane_coord_idx,plane_coord_val,outfpath):
  """Do all the steps for a plane extraction.

  Arguments:

    - infpath = path to the input dx file
    - plane_coord_idx: index of the cartesian axis perpendicular to the plane
        (that is, the cartesian value that is constant in the plane).
        0 for X-plane, 1 for Y-plane, 2 for Z-plane.
    - plane_coord_val: value for the ordinate of the plane, as float
    - outfpath = path to the output pickle file
  
  No return value."""
  in_grid,in_edges=load_dx_convert(infpath)
  out_grid,out_edges=extract_plane(in_grid,in_edges,plane_coord_idx,plane_coord_val)
  save_pickle(out_grid,out_edges,outfpath)
  return
