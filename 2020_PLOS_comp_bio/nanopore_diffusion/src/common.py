"""Classes and functions used in many of the other modules"""

#Standard library
import argparse
from functools import reduce
from itertools import chain
import operator
import os.path as osp
import pickle

#Site packages
import yaml
##from ruamel.yaml import YAML
##yaml=YAML()
try:
  from doit.tools import config_changed
except ImportError:
  #dummy version of config_changed to be used when doit is not available
  def config_changed(arg):
    return arg

#Constants
##pickle_protocol = 4 #The newest protocol, requires python 3.4 or above.
pickle_protocol = 2 #For compatibility with the ancient Python 2

#-------------------------------------------------------------------------------
#Get a value from a nested dictionary

def nested_location(obj,loc):
  """For nested dictionary obj, get item at location specified by the sequence loc

  >>> obj={'a':{'b':{'c':11,'d':99}}}
  >>> loc=['a','b','c']
  >>> nested_location(obj,loc)
  11
  """
  loclist=list(loc)
  return reduce(operator.getitem, [obj]+loclist)

def nested_to_str(obj):
  """For nested dictionary obj, return a consistent string"""
  #Get pairs in correct order, with keys converted to strings
  pairs1=[(repr(k),v) for k,v in obj.items()]
  pairs1.sort(key=lambda t: t[0])
  #Get list of sorted key strings alone
  skeys=[p[0] for p in pairs1]
  #Convert values to strings while maintaining order
  svals=[]
  for k,v in pairs1:
    if hasattr(v,'items'):
      svals.append(nested_to_str(v))
    else:
      svals.append(repr(v))
  items_list=[k+': '+v for k,v in zip(skeys,svals)]
  return '{%s}'%', '.join(items_list)

#And from a nested set of attributes ("drilling down")
def nested_attributes(obj,loc):
  """Liked nested_location, but with attributes instead of dictionary keys"""
  loclist=list(loc)
  return reduce(getattr, [obj]+loclist)

#-------------------------------------------------------------------------------
#Standard work with yaml files

def readyaml(fpath):
  "read object from yaml file"
  with open(fpath,'r') as fp:
    dat=fp.read()
    obj=yaml.load(dat)
  return obj

def readyaml_multidoc(fpath):
  "Return multiple documents from a yaml file as a list"
  with open(fpath,'r') as fp:
    dat=fp.read()
    gen=yaml.load_all(dat)
  return [obj for obj in gen]

def readyaml_gen(fpath):
  "return a generator for the contents of the yaml file"
  with open(fpath,'r') as fp:
    dat=fp.read()
  gen=yaml.load_all(dat)
  return gen

def writeyaml(obj,fpath):
  "write object to yaml file, overwriting"
  with open(fpath,'w') as fp:
    yaml.dump(obj,fp)
  return

#-------------------------------------------------------------------------------
#Standard work with pickle files

def readpickle(fpath):
  "read object from pickle file"
  with open(fpath, 'rb') as fp:
    obj=pickle.load(fp)
  return obj

def writepickle(obj,fpath):
  "write object to pickle file, overwriting"
  with open(fpath, 'wb') as fp:
    pickle.dump(obj,fp,pickle_protocol)
  return

#-------------------------------------------------------------------------------
#ParameterSet

def getbasename(fpath):
  "Return the stem name of the specified file (no directory, no extension)"
  return osp.splitext(osp.basename(fpath))[0]

def add_file_info(d,fpath):
  """Add info from the file path to the dictionary"""
  defaults={'sourcefile':fpath, 'basename':getbasename(fpath)}
  for k,v in defaults.items():
    if k not in d.keys():
      d[k]=v

class ParameterSet(object):
  """A base class for defining sets of related parameters.
  Each subclass should use __slots__ to define its parameters.
  This is not intended as a method for storing complicated objects;
  all attributes should have values that are numbers, strings, sequences, or dictionaries
  whose items follow the same rules.
  
  Attributes:
  
    - basename = a base name which may be shared by a group of instances (usually the stem name of the source file)
    - sourcefile = file from which the object was read, if any

  Frequently, derived classes will add slots for some of the following attributes, which this class supports:

    - _inputfile_attrs: list of attributes containing input file paths
    - _more_inputfiles: list of additional input file paths not contained in other attributes
    - _outputfile_attrs: list of attributes containing output file paths
    - _more_outputfiles: list of additional output file paths not contained in other attributes
    - _folders: dictionary {attribute name: folder path} used to find the folders for prepending to file names specified in other attributes
    - _required_attrs: list of attribute names that must be defined when the object is first loaded
    - _config_attrs: list of attribute names that contain the "configuration" of the object, to be included when doing a check for changes to configuration
    - _taskname_src_attr: string that is the name of the attribute that contains the task name for the object
    - _child_attrs: list of attributes that contain other subclasses of ParameterSet,
      which should be visted recursively when collecting all the input and output files"""
  __slots__=('basename','sourcefile') #Needed even if empty: without this, a __dict__ object will be created even though subclasses use __slots__
  def __init__(self,**kwd):
    ##self.__dict__.update(kwd) #Using __slots__ means there is no __dict__
    #Load the attributes specified
    for k,v in kwd.items():
      setattr(self,k,v)
    #Check for required attributes that are missing
    if hasattr(self,'_required_attrs'):
      missing = [a for a in self._required_attrs if not hasattr(self,a)]
      assert len(missing)==0, "%s missing required attributes: %s"%(type(self),missing)
  @classmethod
  def from_yaml(cls,fpath):
    """Read ParameterSet from a yaml file.

    Arguments:

      - fpath = path to the yaml file to read in
      
        This must be a single-document yaml file,
        and it is assumed to be structured as a single dictionary at the top level.

    Returns:
    
      - pset = a ParameterSet object as defined by the contents of the yaml file"""
    d=readyaml(fpath)
    add_file_info(d,fpath)
    return cls.from_dict(d)
  @classmethod
  def all_from_yaml(cls,fpath):
    """Generator to read a series of ParameterSet objects from a yaml file.

    Arguments:

      - fpath = path to the yaml file to read in

        This must be a multi-document yaml file,
        and each document is assumed to be structured as a single dictionary at the top level.

    Each call yields:

      - pset = a ParameterSet object as defined by one document in the yaml file"""
    with open(fpath,'r') as fp:
      dat=fp.read()
      gen=yaml.load_all(dat)
    for d in gen:
      add_file_info(d,fpath)
      yield cls.from_dict(d)
  @classmethod
  def select_from_yaml(cls,fpath,criteria):
    """Generator to read ParameterSet objects matching criteria from a yaml file.

    Arguments:

      - fpath = path to the yaml file to read in

        This must be a multi-document yaml file,
        and each document is assumed to be structured as a single dictionary at the top level.

      - criteria = dictionary of criteria

        Each key in this dictionary must have a matching attribute in the ParameterSet.
        Only objects whose attribute values match those specified by this dictionary are yielded.

    Each call yields:

      - pset = a ParameterSet object as defined by one document in the yaml file"""
    for pset in cls.all_from_yaml(fpath):
      for k,v in criteria.items():
        if getattr(pset,k,None) == v:
          yield pset
  @classmethod
  def select_one_from_yaml(cls,fpath,criteria):
    """To return a single ParameterSet object matching criteria from a yaml file

    Arguments:

      - fpath = path to the yaml file to read in

        This must be a multi-document yaml file,
        and each document is assumed to be structured as a single dictionary at the top level.

      - criteria = dictionary of criteria

        Each key in this dictionary must have a matching attribute in the ParameterSet.
        Only objects whose attribute values match those specified by this dictionary are yielded.

    Each call yields:

      - pset = a ParameterSet object as defined by one document in the yaml file"""
    psets=[p for p in cls.select_from_yaml(fpath,criteria)]
    assert len(psets)==1, "Expected only one result, but instead got %d"%len(psets)
    return psets[0]
  def to_yaml(self,fpath):
    """Write ParameterSet to a yaml file.

    Arguments:

      - fpath = path to the yaml file to write

        This will be a single-document yaml file, containing a dictionary (potentially of other dictionaries).
        The file will be overwritten if it already exists.

    No return value."""
    #TODO: consider writing docstring as comments in yaml file
    writeyaml(self.to_dict(),fpath)
    return
  @classmethod
  def from_pickle(cls,fpath):
    """Read ParameterSet from a pickle file.

    Arguments:

      - fpath = path to the pickle file to read in

        This file should be a mapping type at the top level.

    Returns:

      - pset = a ParameterSet object as defined by the contents of the pickle file."""
    d=readpickle(fpath)
    add_file_info(d,fpath)
    return cls.from_dict(d)
  def to_pickle(self,fpath):
    """Write ParameterSet to a pickle file.

    Arguments:

      - fpath = path to the pickle file to write
        The file will be overwritten if it already exists.

    No return value."""
    writepickle(self.to_dict(),fpath)
    return
  @classmethod
  def from_dict(cls,d):
    """Load the object from a dictionary"""
    return cls(**d)
  def _all_slots_iter(self):
    """Return an iterator over all the available slots"""
    return chain.from_iterable(getattr(cls, '__slots__', []) for cls in type(self).__mro__)
  def to_dict(self):
    """Return a dictionary with all the object's attributes.

    Note that changes to this dictionary will not affect the object.

    No arguments.

    Returns the dictionary."""
    return dict([(k,getattr(self,k,None)) for k in self._all_slots_iter()])
  def __repr__(self):
    return nested_to_str(self.to_dict())
  @classmethod
  def from_Namespace(cls,ns):
    return cls(**ns.__dict__)
  def to_Namespace(self):
    """Return an argparse.Namespace object with all the object's attributes.

    Note that changes to the Namespace will not affect the object.

    No arguments.

    Returns the Namespace."""
    return argparse.Namespace(**self.to_dict())
  ##TODO: read and write from ini file
  #Methods supporting task definition for doit
  @property
  def task_definition(self):
    """A doit task definition dictionary appropriate for processing this object."""
    return {'name': self.taskname,
     'file_dep': self.inputfiles,
     'uptodate': [config_changed(self.config)],
     'targets': self.outputfiles,
     'actions': [(self.run,)]}
  @property
  def config_dict(self):
    d=self.to_dict()
    cd=dict([(k,v) for k,v in d.items() if k in self._config_attrs])
    #Don't add configuration of children: look at the output files they generate instead
    #(otherwise changes will cascade even if output files are unaltered)
    ##for childattr in getattr(self,'_child_attrs',[]):
    ##  cd[childattr]=getattr(self,childattr).config
    return cd
  @property
  def config(self):
    """A string representing the configuration of the object, suitable for use by doit.tools.config_changed."""
    # return(str(self.config_dict))
    return(yaml.dump(self.config_dict))
  @property
  def taskname(self):
    """A string representing the task name in doit"""
    return getattr(self,self._taskname_src_attr)
  @property
  def inputfiles(self):
    """A list of all the inputfiles related to this object"""
    ifl=self.file_list('_inputfile_attrs')
    ifl+=getattr(self,'_more_inputfiles',[])
    for childattr in getattr(self,'_child_attrs',[]):
      ifl += getattr(self,childattr).inputfiles
    return ifl
  @property
  def outputfiles(self):
    """A list of all the outputfiles related to this object"""
    ofl = self.file_list('_outputfile_attrs')
    ofl+=getattr(self,'_more_outputfiles',[])
    for childattr in getattr(self,'_child_attrs',[]):
      ofl += getattr(self,childattr).outputfiles
    return ofl
  def file_list(self,top_attr):
    """Return the files whose relative paths are provided in top_attr"""
    return [self.full_path(attr) for attr in getattr(self,top_attr,[])]
  def full_path(self,attr):
    """Return the full path string associated with a specific attribute"""
    return osp.join(getattr(self,'_folders',{}).get(attr,''),getattr(self,attr))

#-------------------------------------------------------------------------------
#Common command-line usage

def run_cmd_line(program_description,input_file_description,objtype,process_method="run",other_selection=None,f_args=None):
  """Perform common command line processing.

  Many of the modules use a similar process for running from the command line:

  - read an multi-document input yaml file
  - load each document into an object of some type
  - select specific objects using the "--select" command line argument, and/or other means
  - run a function on each selected object

  This function implements that common process.

  Arguments:

    - program_description = string containing the help for the program to run
    - input_file_desscription = string containing help for the input file to process
    - objtype = type all documents in the input file should be loaded into

      This is assumed to be a subclass of ParameterSet.
      (At a minimum, it must have an all_from_yaml method)

    - process_method = optional string, name of object's method to be called
    - other_selection = optional dictionary with additional requirements for objects to process

      {attribute name: sequence of allowed values}

    - f_args = optional list of positional arguments for calling the process_method

  No return value."""

  #Parse arguments
  parser = argparse.ArgumentParser(description=program_description)
  parser.add_argument('params_yaml', help=input_file_description)
  parser.add_argument('--select',nargs="+",help="""Only process selected elements of the input yaml file.
    This option must be followed by an attribute name, and then a sequence of values for that attribute.
    Only those entries in the yaml file where the attribute matches one of these values will be processed.
    If the attribute in question will be cast to a string, for comparison to the string read from the command line.
    Note, in particular, this means you may need to explicitly provide escaped quotes in some cases.""")
  cmdline=parser.parse_args()
  
  #Set up dictionary for object selection by attribute
  if other_selection is None:
    requirements = {}
  else:
    requirements=dict(other_selection)
  if cmdline.select is not None:
    requirements[cmdline.select[0]]=cmdline.select[1:]

  #Go
  common_run(cmdline.params_yaml,objtype,process_method,requirements,f_args)


def common_run(input_file,objtype,process_method="run",requirements=None,f_args=None):

  #Check that input file exists
  assert osp.isfile(input_file), "Parameter definition file does not exist: %s"%input_file

  #Dictionary for object selection by attribute
  if requirements is None:
    requirements={}

  #Any other positional arguments?
  if f_args is None:
    f_args = []

  #Read in the yaml file
  allobjs=objtype.all_from_yaml(input_file)

  #Process documents
  for obj in allobjs:
    #Is this document selected? (if no selection list provided, process all documents)
    selected=True
    for selattr,allowed in requirements.items():
       selected = selected and (str(getattr(obj,selattr,None)) in allowed)
       if not selected:
         break
    if selected:
      getattr(obj,process_method)(*f_args)
