"""Functions useful in converting xvg files to pandas dataframes

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os.path as osp

#Site packages
import numpy as np
import pandas as pd

#Constants
IGNORECHARS='#@'
SET_SEPCHAR='&'
traj_headings=['t','x','y','z']
traj_structure=[4]
acf_headings=[['t','x'], ['t','y'], ['t','z']]
acf_structure=[2,2,2]

class FoundNewSet(Exception):
  pass

def provide_sets(infpath):
  """Generator yielding all the data sets in the given xvg file.

  At the start of a new data set (other than the first one),
  an instance of ``FoundNewSet`` is returned.

  Arguments:

    - infpath = path to input file, as string
  
  Yields:

    - entries = each row of each data set in the file"""
  #Open the file, and close properly when done
  with open(infpath,'r') as fp:
    #Read each line separately
    for line in fp.readlines():
      #Ignore lines starting with specified characters:
      if line[0] not in IGNORECHARS:
        #Raise FoundNewSet if a new set starts here
        if line[0]==SET_SEPCHAR:
          yield FoundNewSet()
        #Otherwise, yield another row of data
        else:
          yield line.split()
  #Done
  return

##TODO: a version of this function that works on output of read_all_sets, instead of reading the file itself
def get_structure(infpath):
  """Get the structure of the data sets in the given xvg file.

  Arguments:

    - infpath = path to input file, as string
  
  Returns:

    - setcounts = list with one number for each set in the file:
                  the number of entries in each row for that set
  
  If a set is found to have different numbers of entries in different rows,
  an AssertionError is raised."""
  #Initializations
  setcounts=[]
  thisset=None
  #Go through the sets and rows
  for entries in provide_sets(infpath):
    if isinstance(entries,FoundNewSet):
      setcounts.append(thisset)
      thisset=None
    else:
      if thisset is None:
        thisset = len(entries)
      else:
        assert thisset == len(entries), "Previous rows had %d entries per row, this one has %d: %s"%(thisset, len(entries),str(entries))
  #Done
  if thisset is not None:
    setcounts.append(thisset)
  return setcounts

def read_all_sets(infpath,typefunc=float):
  """Read all the data sets in the given xvg file.

  Arguments:

    - infpath = path to input file, as string
    - typefunc = optional, type/class to use for each value, defaults to ``float``
  
  Returns:

    - allsets = list of data sets in the file, each set a list of rows, each row a list of values"""
  #Initializations
  allsets=[]
  setnow=[]
  #Go through the sets and rows
  for entries in provide_sets(infpath):
    if isinstance(entries,FoundNewSet):
      allsets.append(setnow)
      setnow=[]
    else:
      vals=[typefunc(itm) for itm in entries]
      setnow.append(vals)
  #Done
  if len(setnow)>0:
    allsets.append(setnow)
  return allsets

def sets_to_frames(allsets,columns):
  """Convert a sequence of data sets to a sequence of frames

  Arguments:

    - allsets = list of data sets, each data set a list of rows of values
    - columns = a sequence with one entry for each data set: a sequence of column headings

  The structure of ``columns`` must match the structure of ``allsets``.

  Returns:

    - allframes = a list of pandas dataframes"""
  #Error checking setup
  outer_err_template = "Must have one group of column headings for each data set: got %d data sets, %d heading groups"
  inner_err_template = "Data Set index %d: Must have one column heading for each entry in the rows: %d entries, %d column headings"
  #Check that outer structure matches
  outer_err_tuple = (len(allsets),len(columns))
  assert outer_err_tuple[0]==outer_err_tuple[1], outer_err_template%outer_err_tuple
  #Initialization
  allframes=[]
  #Go through each set
  for setnum,thisset in enumerate(allsets):
    thisheader=columns[setnum]
    #Check that inner structure matches
    inner_err_tuple = (setnum,len(thisset[0]),len(thisheader))
    assert inner_err_tuple[1]==inner_err_tuple[2], inner_err_template%inner_err_tuple
    #Create DataFrame
    df=pd.DataFrame(thisset,columns=thisheader)
    allframes.append(df)
  #Done
  return allframes

def combine_frames(framelist,matching):
  """Combine DataFrames that have matching columns

  Arguments:

    - framelist = sequence of DataFrames to be combined
    - matching = list of column labels where data should match for all rows
  
  Returns:

    - outframe = combined DataFrame

  Only one copy of the matching columns are included in the output DataFrame."""
  #Initial checks
  assert len(framelist)>1, "You gave me one DataFrame. What did you want me to combine it with?"
  assert len(matching)>0, "Must have at least one column matching between DataFrames"
  #Copy the first frame
  outframe = framelist[0].copy()
  #Go through the other frames
  for otherframe in framelist[1:]:
    #Confirm that the matching columns do match
    for mcol in matching:
      match_result = all(outframe.loc[:,mcol]==otherframe.loc[:,mcol])
      assert match_result, "Column %s does not match for frame indices 0 and %d"%(str(mcol),frameno)
    #Get a list of the columns to be added
    addcols=[c for c in otherframe.columns if c not in matching]
    #Add those columns
    for colname in addcols:
      outframe[colname] = otherframe[colname]
  #Done
  return outframe

def load_acf_xvg(fpath):
  """Load an xvg file containing an ACF"""
  assert osp.isfile(fpath)
  corrsets=read_all_sets(fpath)
  ##TODO: check structure of the file
  corrframes=sets_to_frames(corrsets,acf_headings)
  corrdf=combine_frames(corrframes,['t'])
  return corrdf

def load_traj_xvg(fpath):
  """Load an xvg file containing a trajectory"""
  assert osp.isfile(fpath)
  trajsets=read_all_sets(fpath)
  ##TODO: check structure of the file
  trajframes=sets_to_frames(trajsets,[traj_headings])
  return trajframes[0]
