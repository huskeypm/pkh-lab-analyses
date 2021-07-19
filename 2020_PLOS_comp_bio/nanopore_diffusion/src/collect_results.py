"""For collecting results from multiple simulations"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp

#Site packages
import pandas as pd

#Local
import folderstructure as FS
import common

#Constants
collected_df_fname='collected_results.pkl.gz'

def get_columns(d,exclusions):
  """Return list of all keys where the values are integers, floats, or strings,
  and call recursively on any value that has its own 'items' attribute.
  Sequences are now accepted as well, with an index added to the name of their parent"""
  cols=[]
  newcols=[]
  for k,v in d.items():
    if not k in exclusions:
      if type(v)==int or type(v)==float or type(v)==str:
        if not k in cols:
          cols.append(k)
      elif hasattr(v,'index'):
        for idx,itm in enumerate(v):
          sub_dict={k+"_%d"%idx: itm}
          newcols+=get_columns(sub_dict,exclusions)
      elif hasattr(v,'items'):
        newcols=get_columns(v,exclusions)
      for c in newcols:
        if c not in cols:
          cols.append(c)
      newcols=[]
  return cols

def get_all_columns(dlist,exclusions):
  """Return the superset of columns for each d in dlist"""
  columns=get_columns(dlist[0],exclusions)
  for d in dlist[1:]:
    newcols=get_columns(d,exclusions)
    for c in newcols:
      if not c in columns:
        columns.append(c)
  return columns

def flatdict(d,cols,name,exclusions):
  """Flatten a potentially nested dictionary so it can be added to a DataFrame with the specified columns."""
  fd={}
  subdict={}
  for k,v in d.items():
    if not k in exclusions:
      if k in cols and (type(v)==int or type(v)==float or type(v)==str):
        fd[k]=v
      elif hasattr(v,'index'):
        subdict={}
        for idx,itm in enumerate(v):
          newk=k+"_%d"%idx
          newname=name+'->'+newk
          sd={newk: itm}
          subdict.update(flatdict(sd,cols,newname,exclusions))
      elif hasattr(v,'items'):
        newname=name+'->'+k
        subdict=flatdict(v,cols,newname,exclusions)
      for newk, newv in subdict.items():
        if newk in fd.keys():
          assert newv==fd[newk], "In %s, unequal assignments to %s: previously had %s, but %s wants to change it to %s"%(name,str(newk),str(fd[newk]),str(k),str(newv))
        else:
          fd[newk]=newv
      subdict={}
  return fd

def dicts_to_dataframe(alldicts,exclusions):
  """Create a pandas dataframe from an iterable of dictionaries.

  Arguments:

    - alldicts = iterable of dictionaries
    - exclusions = list of keys to exlcude

  Return value:

    - df = pandas dataframe

  For each dictionary:

    - Anything that is a number or string is added directly.
    - Anything that is a dictionary has its items treated the same way.
      (More specifically, anything that has an 'items' attribute.)
    - Everything else is ignored, including any sequences."""
  #Get the list of columns for the dataframe, and make it works for all the dictionaries
  columns=get_all_columns(alldicts,exclusions)
  #Set up dataframe with required columns
  df=pd.DataFrame(columns=columns)
  #Add the data to the dataframe
  for d in alldicts:
    fd=flatdict(d,columns,FS.infofile,exclusions)
    df=df.append(fd,ignore_index=True)
  return df
  
class ResultsCollector(common.ParameterSet):
  __slots__=('modellist','exclusions','input_files','outfpath')
  _required_attrs=['basename','modellist']
  _config_attrs=['modellist','exclusions']
  _outputfile_attrs=['outfpath']
  _taskname_src_attr=['taskname']
  
  def __init__(self,**kwd):
    #Initialization from base class
    super(ResultsCollector, self).__init__(**kwd)
    #Get the output file
    self.outfpath=osp.join(FS.postprocfolder,self.basename,collected_df_fname)
    #Get the input files
    basedir=osp.join(FS.solnfolder,self.basename)
    self.input_files=[]
    for modelparams in self.modellist:
      if modelparams.basename==self.basename:
        #Add this one to the list
        self.input_files.append(osp.join(basedir,modelparams.modelname,FS.infofile))

  @property
  def taskname(self):
    return self.basename+":collection"
  
  @property
  def _more_inputfiles(self):
    return self.input_files

  def run(self):
    print(self.outfpath)
    #Make sure the parent folder exists
    if not osp.isdir(osp.split(self.outfpath)[0]):
      os.makedirs(osp.split(self.outfpath)[0])
    #Read in all the info files
    alldicts = [common.readyaml(fp) for fp in self.input_files]
    #Convert to dataframe and save
    df = dicts_to_dataframe(alldicts,self.exclusions)
    df.to_pickle(self.outfpath)
    