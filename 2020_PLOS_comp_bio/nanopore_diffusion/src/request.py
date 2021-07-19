"""Define the base class for all Requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
from itertools import chain
import pickle

#Site packages
import yaml
##from ruamel.yaml import YAML
##yaml=YAML()

#Local
import filepath

#Constants
##pickle_protocol = 4 #The newest protocol, requires python 3.4 or above.
pickle_protocol = 2 #For compatibility with the ancient Python 2

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

class ParameterSet(object):
  """A base class for defining sets of related parameters.
  Each subclass should use __slots__ to define its parameters.
  This is not intended as a method for storing complicated objects;
  all attributes should have values that are numbers, strings, sequences, or dictionaries
  whose items follow the same rules.
  
  Basically, this is just a namespace that knows how to read and write itself to/from file(s) of different types.
  
  Attributes:
  
    - sourcefile = file from which the object was read, if any

  """
  __slots__=('__sourcefile__',) #Needed even if empty: without this, a __dict__ object will be created even though subclasses use __slots__
  def __init__(self,**kwargs):
    ##self.__dict__.update(kwd) #Using __slots__ means there is no __dict__
    #Load the attributes specified
    for k,v in kwargs.items():
      setattr(self,k,v)
    #Define the sourcefile if not already defined
    if '__sourcefile__' not in kwargs.keys():
      self.__sourcefile__ = None
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
    return cls.from_dict(**ns.__dict__)
  def to_Namespace(self):
    """Return an argparse.Namespace object with all the object's attributes.

    Note that changes to the Namespace will not affect the object.

    No arguments.

    Returns the Namespace."""
    return argparse.Namespace(**self.to_dict())
  ##TODO: read and write from ini file
  @classmethod
  def all_from_yaml(cls,fpath):
    """Generator to read a series of objects from a yaml file.

    Arguments:

      - fpath = path to the yaml file to read in, as a pathlike object

        This must be a multi-document yaml file,
        and each document is assumed to be structured as a single dictionary at the top level.

    Each call yields:

      - pset = a ParameterSet object as defined by one document in the yaml file"""
    with open(fpath,'r') as fp:
      dat=fp.read()
      gen=yaml.load_all(dat)
    for d in gen:
      d.update({'__sourcefile__':fpath})
      yield cls.from_dict(d)
  @classmethod
  def one_from_yaml(cls,fpath):
    """Read only a single instance of objects from a yaml file.

    Arguments:

      - fpath = path to the yaml file to read in
      
        This must be a single-document yaml file,
        and it is assumed to be structured as a single dictionary at the top level.

    Returns:
    
      - pset = a ParameterSet object as defined by the contents of the yaml file"""
    for pset in cls.all_from_yaml(fpath):
      return pset #ugly, but effective
  def to_yaml(self,fpath):
    """Write a single instance to a yaml file.

    Arguments:

      - fpath = path to the yaml file to write, as a pathlike objeect

        The method will prepend a document separator to the object data,
        and will "append" to an existing file, allowing multiple documents to be written to the same file.

    No return value."""
    #TODO: consider writing docstring as comments in yaml file
    obj=self.to_dict()
    with open(fpath,'a') as fp:
      fp.write("---\n")
      yaml.dump(obj,fp)
    return
  ##TODO: both pickle methods should support sequences of objects as well as single objects
  @classmethod
  def from_pickle(cls,fpath):
    """Read from a pickle file.

    Arguments:

      - fpath = path to the pickle file to read in

        This file should be a mapping type at the top level.

    Returns:

      - pset = a ParameterSet object as defined by the contents of the pickle file."""
    with open(fpath, 'rb') as fp:
      d=pickle.load(fp)
    d.update({'__sourcefile__':fpath})
    return cls.from_dict(d)
  def to_pickle(self,fpath):
    """Write to a pickle file.
    
    This method doesn't prevent the ParameterSet from having entries which are instances of a class,
    but this is highly discouraged.
    Entries should generally consist only of the same data types supported by YAML.

    Arguments:

      - fpath = path to the pickle file to write
        The file will be overwritten if it already exists.

    No return value."""
    obj=self.to_dict()
    with open(fpath, 'wb') as fp:
      pickle.dump(obj,fp,pickle_protocol)
    return
    
class Request(ParameterSet):
  __slots__=('requestclass')
  ##TODO: everything
