"""Get and set values for not just simple attributes or keys, but dotted paths of attributes/keys"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
from copy import deepcopy

#Site packages
import ruamel.yaml

#This package
from . import yaml_manager

#Helper functions

def to_sequence(dpath):
  if isinstance(dpath,str):
    seq=dpath.split('.')
  else:
    seq=dpath
  return seq

def drill_down(obj,seq):
  head=seq[:-1]
  tail=seq[-1]
  if hasattr(obj,'get_nested'):
    parent=obj.get_nested(head)
  else:
    parent=get_nested(obj,head)
  return parent,tail

def get_nested(obj,dpath):
  """Return the value from the specified attribute/key/index path
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
        The path may contain attributes and dictionary keys, with no need to distinguish between them.
        List indices are also allowed, but only for sequence arguments,
        as strings are not cast to other data types.
  
  Returns the requested data."""
  nxt=obj
  seq=to_sequence(dpath)
  for name in seq:
    if isinstance(name,str) and hasattr(nxt,name):
      nxt = getattr(nxt,name)
    else:
      try:
        #Note that this will work for both lists and dictionaries
        nxt=nxt.__getitem__(name)
      except Exception as the_err:
        raise KeyError('Invalid path %s: No attribute, key, or index %s'%(dpath,name)) from the_err
  return nxt

def get_nested_default(obj,dpath,default=None):
  """Just like get_nested, but return a default object if any part of the path is invalid"""
  try:
    return get_nested(obj,dpath)
  except:
    return default

def set_nested(obj,dpath,val):
  """Set the value at the specified attribute/key/index path
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
    - val = value to assign
  
  No return value."""
  seq=to_sequence(dpath)
  parent,tail=drill_down(obj,seq)
  if hasattr(parent,'__setitem__'):
    parent.__setitem__(tail,val)
  else:
    setattr(parent,tail,val)
  return

def update_nested(obj,dpath,val):
  """Call the ``update`` method of the object at the specified attribute/key/index path

  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
    - val = value to be used as the argument to ``update``

  No return value."""
  seq=to_sequence(dpath)
  if hasattr(obj,'get_nested'):
    target=obj.get_nested(seq)
  else:
    target=get_nested(obj,seq)
  target.update(val)

def new_odict(obj,dpath):
  """Create a new OrderedDict instance suitable for later use with set_nested
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
  
  No return value."""
  seq=to_sequence(dpath)
  parent,tail=drill_down(obj,seq)
  setattr(parent,tail,odict())
  return

#Classes

class WithNested(dict):
  """A very basic class that loads and sets nested attributes, and supports reading and writing itself to/from yaml.

  This class is derived from ``dict``, and unifies access to data through dictionary and object methods.
  That is, its attributes and keys are identical.

  It also has the ability to track all members of its class and subclasses with a name in a global dictionary.
  
  Supported attributes:
  
    - name: name for the object, as a string
    - store_globally: boolean, True to add to the global dictionary of named objects
  
  Note that if no name is provided, the object cannot be stored globally.
  
  Example uses:
  
  >>> bunch=WithNested(a=1,b=2,c=3)
  >>> bunch
  {'a': 1, 'b': 2, 'c': 3}
  >>> bunch.a
  1
  >>> bunch.b=-99
  >>> bunch['b']
  -99"""
  __allnames__=odict()
  def __init__(self,**kwargs):
    #Unify access through attributes and keys
    self.__dict__=self
    #Load the attributes specified
    for k,v in kwargs.items():
      setattr(self,k,v)
    #If requested, add this object to the global dictionary
    if getattr(self,'store_globally',False):
      assert getattr(self,'name',None) is not None, "Cannot add unnamed object to global names dictionary"
      assert self.name not in self.__allnames__.keys(), "Attempted to overwite existing global name: %s"%self.name
      self.__allnames__[self.name]=self
  def to_dict(self,recursive=True):
    """Return a dictionary with all the object's attributes.

    Note that changes to this dictionary will not affect the object.

    No arguments.

    Returns the dictionary."""
    d={}
    for attr,itm in self.items():
      if recursive and hasattr(itm,'to_dict') and callable(itm.to_dict):
        d[attr]=itm.to_dict()
      else:
        d[attr]=itm
    return d
  def __getstate__(self):
    """Used for pickling, and converting to yaml"""
    return self.to_dict(recursive=True)
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)
  # def get(self,loc,default=None):
  #   if default is None:
  #     return self.__dict__.get(loc)
  #   else:
  #     return self.__dict__.get(loc,default)
  def get_copy(self):
    """Return a deep copy of this instance"""
    return deepcopy(self)
  def get_nested(self,dpath):
    """Return the value from the specified attribute/key/index path"""
    return get_nested(self,dpath)
  def get_nested_default(self,dpath,default=None):
    """Return the value from the specified attribute/key/index path, with default value if not present"""
    return get_nested_default(self,dpath,default)
  def set_nested(self,dpath,val):
    """Set the value at the specified attribute/key/index path"""
    set_nested(self,dpath,val)
    return
  def update_nested(self,dpath,val):
    """Call ``update`` on the object at the specified attribute/key/index path"""
    update_nested(self,dpath,val)
    return
  def new_odict(self,dpath):
    """Set up a new OrderedDict for later use with set_nested"""
    new_odict(self,dpath)
    return
  def get_stored(self,candidate):
    """Get the value from the specified storage location

    This is used to process instances of the ``Stored`` class.
    
    If the candidate is an instance of ``Stored``, get what it points to.
    Otherwise, return the candidate unaltered."""
    if isinstance(candidate, Stored):
      return self.get_nested(candidate.attrloc)
    else:
      return candidate

class Stored(object):
  """Provide the location of a value stored elsewhere in memory

    (See ``WithNested.get_stored`` for internal usage.)
  
    Attributes:
    
      - attrloc = basically a dpath argument for use with get_nested and set_nested
  """
  def __init__(self,attrloc):
    self.attrloc=attrloc
  @classmethod
  def from_yaml(cls, constructor, node):
    if type(node)==ruamel.yaml.nodes.SequenceNode:
      return cls([constructor.construct_scalar(itm) for itm in node.value])
    else:
      #type(node)==ruamel.yaml.nodes.ScalarNode
      return cls(node.value)
  @classmethod
  def to_yaml(cls, representer, node):
    return representer.represent_scalar("!"+cls.__name__,str(node.attrloc))

#Register for loading from yaml
yaml_manager.register_classes([Stored])
