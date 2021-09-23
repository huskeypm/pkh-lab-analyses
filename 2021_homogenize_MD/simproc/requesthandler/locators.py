"""Define data locators

Locators are objects that can return a path to a particular file.
First, an overall file structure is defined.
Within this file structure, individual paths are computed
based on the parameters of the request that needs them.

This module has a variable ``DATAFOLDER`` which specifies the path of the folder to work in.
All locators return paths that start at this folder.

From within python, you can simply set ``DATAFOLDER`` as a variable,
and all locators will return paths based on the new location.

For scripts, use the environment variable ``DATAFOLDER`` to tell the program where this folder is.

If you want to see the current list of locators, just look at ``folder_structure``."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
import os

#Site packages

#This package
from . import filepath
from . import yaml_manager

#Constants
datafolder_environ='DATAFOLDER'

#Get path to the top folder
if datafolder_environ in os.environ.keys():
  # DATAFOLDER=Path(osp.normpath(osp.abspath(os.environ[datafolder_environ])))
  DATAFOLDER=filepath.Path(os.environ[datafolder_environ]).expanduser().resolve()
else:
  modpath=filepath.Path(__file__)
  srcfolder=modpath.parent.parent
  DATAFOLDER=srcfolder.parent / 'data'

class DataFile(object):
  """File location relative to ``DATAFOLDER``
  
  The location of ``DATAFOLDER`` is specified by the environment variable ``DATAFOLDER``.
  If this environment variable does not exist, a default value is provided."""
  def __init__(self,*args,**kwargs):
    self.subpath=filepath.Path(*args,**kwargs)
  def path(self,reqname):
    return DATAFOLDER / self.subpath
  @classmethod
  def from_yaml(cls, constructor, node):
    return cls(node.value)
  @classmethod
  def to_yaml(cls, representer, node):
    return representer.represent_scalar("!"+cls.__name__,str(node.subpath))
  def __repr__(self):
    return self.__class__.__name__+"("+repr(self.subpath)+")"

class NameDelegator(object):
  """Act like a locator, but with an alternative request name."""
  def __init__(self,req,loc):
    self.req=req
    self.loc=loc
  def path(self,wrong_reqname):
    return self.loc.path(self.req)

class Delegator(object):
  """Act like a locator, but delegate to another request"""
  def __init__(self,reqattr,loc):
    self.reqattr=reqattr
    self.loc=loc
  def render(self,parent):
    req=parent.get_nested(self.reqattr)
    return req.render(self.loc)

def locator_factory(ltype):
  """Factory function to return a locator class for a given name
  
  The locator class returned will construct file paths as follows:
  
    - all paths begin with ``locators.DATAFOLDER``
    - then, one subdirectory is added for each element of the folder structure definition sequence for the locator (see below)
    - finally, the subpath used to initialize the locator is appended
  
  The folder structure definition sequence defines subdirectories as follows:
  
    - any element which is a string adds the string as a folder name
    - any element which is an integer adds the matching element in the dotted name of the request
    - integers beyond the last element of the dotted name are simply ignored.
  
  Arguments:
  
    - ltype = name of locator class, as string
  
  Returns:
  
    - lclass = locator class"""
  class lclass(DataFile):
    """File location as specified by the name of its parent request and the current folder structure settings"""
    def __init__(self,*args,**kwargs):
      self.subpath=filepath.Path(*args,**kwargs)
    def path(self,reqname):
      namelist=reqname.split('.')
      specifier=folder_structure[ltype]
      out=DATAFOLDER
      for itm in specifier:
        if type(itm) is int:
          if itm < len(namelist):
            out /= namelist[itm]
        else:
          out /= itm
      out /= self.subpath
      return out
  lclass.__name__=ltype
  return lclass

class FolderStructure(odict):
  """For storing an expected folder structure.
  
  This dictionary is intended to have the following organization: {locator_class_name: [structure_defintion_sequence]}  
  
    - each key is a string, matching the name of a locator class
    - the structure definition sequence indicates how to create the path for files using the locator class
    
  See the locator class factory for explanation of how the structure definition sequence determines the file path.""" 
  def update(self,**kwargs):
    #Before we actually change the folder structure, we need to check for new locator classes
    #All existing locator types
    existing=self.keys()
    #Newly defined types
    newtypenames=[lt for lt in kwargs.keys() if lt not in existing]
    newtypes=[]
    for ltype in newtypenames:
      #Create the locator
      lclass=locator_factory(ltype)
      globals()[ltype]=lclass
      newtypes.append(lclass)
    #Register the new locators for reading from yaml
    yaml_manager.register_classes(newtypes)
    #Update the expected folder structure
    super(FolderStructure, self).update(**kwargs)

#Folder structure singleton
folder_structure=FolderStructure()

class UpdateFolderStructure(object):
  """Apply changes to the folder structure used to locate files
  
  This isn't really a class. It's a function.
  But when loading from yaml, there isn't a way to call a function directly.
  You can only load classes.
  Hence, this is an object that simply calls another function when initialized,
  and then does nothing else ever.
  
  This is also not a request: its action needs to be taken at initialization,
  not when requests are run, so that the locators in the input can make use of
  the newly defined folder structure.
  Note that this means its attributes are not validated like a request's are."""
  def __init__(self,**kwargs):
    folder_structure.update(**kwargs)
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)
  def __getstate__(self):
    """Used for pickling, and writing to yaml"""
    return folder_structure

# class DumpFolderStructure(object):
#   """Write the current folder structure to a yaml output file

#   This isn't really a class. It's a function.
#   But when loading from yaml, there isn't a way to call a function directly.
#   You can only load classes.
#   Hence, this is an object that simply calls another function when initialized,
#   and then does nothing else ever.
  
#   This is also not a request: its action is taken at initialization,
#   not when requests are run.
#   That's a questionable decision, but the way I decided to do it for now.
#   Note this means its attributes are not validated like a request's are.
  
#   The only input argument is `outfile` which must be a string or Path, NOT a locator,
#   because this isn't a request.
#   The path will be assumed to be relative to DATAFOLDER, as it exists at the time of calling."""
#   def __init__(self,outfile):
#     self.outfile=outfile
#     out_relpath=filepath.Path(outfile)
#     out_abspath=DATAFOLDER / out_relpath
#     out_relpath.assure_dir()
#     d=dict(folder_structure.items())
#     yaml_manager.writefile_flow(d,str(out_abspath))
#   def __setstate__(self,state):
#     """Used for unpickling, and loading from yaml"""
#     self.__init__(**state)
#   def __getstate__(self):
#     """Used for pickling, and writing to yaml"""
#     return {'outfile':self.outfile}

class SetDataFolder(object):
  """Set the DATAFOLDER
  
  As with `UpdateFolderStructure`, this is used as a function.
  It allows DATAFOLDER to be set from within a yaml input file.
  
  It's not a request, for the same reason that `UpdateFolderStructure`
  is not a request either."""
  def __init__(self,**kwargs):
    global DATAFOLDER
    if len(kwargs)>0:
      targpath=filepath.Path(kwargs['datafolder'])
      targpath=targpath.expanduser()
      if not targpath.is_absolute():
        #Relative path, assumed to be relative to the yaml file location (if one is being loaded)
        #If you're not loading a yaml file, it will just use python's current directory
        if len(yaml_manager.now_loading)>0:
          yamlpath_str=yaml_manager.now_loading[-1]
          yamlpath=filepath.Path(yamlpath_str,isfile=True)
          yamldir=yamlpath.folder_path
          targpath=yamldir/targpath
      if kwargs.get('resolve',True):
        targpath=targpath.resolve()
      DATAFOLDER=targpath
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)
  def __getstate__(self):
    """Used for pickling, and writing to yaml"""
    global DATAFOLDER
    return {'datafolder':DATAFOLDER}

#Register for loading from yaml
# yaml_manager.register_classes([filepath.Path,DataFile,NameDelegator,Delegator,UpdateFolderStructure,DumpFolderStructure,SetDataFolder])
yaml_manager.register_classes([filepath.Path,DataFile,NameDelegator,Delegator,UpdateFolderStructure,SetDataFolder])
