"""Command-line support for running requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import importlib
import sys

#Site packages
import yaml
##from ruamel.yaml import YAML
##yaml=YAML()

#Local
import new_folderstructure as FS
import filepath
import request

def obtain_class(requestclass):
  """Return the class specified by a string
  
  The string format is "module_name:class_name".
  """
  modname,classname=requestclass.split(':')
  if modname not in sys.modules.keys():
    importlib.import_module(modname)
  rclass=getattr(sys.modules[modname],classname)
  return rclass

class TotipotentRequest(request.Request):
  """A type of request that can turn itself into any other request defined in a module that can be loaded

  Supported Attributes:
    These are class attributes that Request subclasses may provide, but do not have to.
    
    - _child_attrs = sequence of attribute names, for attributes containing individual child Requests
    - _child_seq_attrs = sequence of attribute names, for attributes containing sequences of child Requests
  """
  @classmethod
  def build_request(cls,**kwargs):
    #Get the class to build
    classname=kwargs['requestclass']
    rclass=obtain_class(classname)
    #Recurse through child requests in the dictionary structure
    for attrname in getattr(rclass,'_child_attrs',[]):
      if attrname in kwargs.keys():
        kwargs[attrname]=cls.build_request(**kwargs[attrname])
    for attrname in getattr(rclass,'_child_seq_attrs',[]):
      if attrname in kwargs.keys():
        kwargs[attrname]=[cls.build_request(**itm) for itm in kwargs[attrname]]
    #Initialize and return object
    return rclass(**kwargs)
  @classmethod
  def from_dict(cls,d):
    """Load the object from a dictionary"""
    return cls.build_request(**d)

class RequestFileRequest(request.Request):
  """Request to run all of the requests listed in the specified files
  
  Attributes:
  
    - requestfiles: sequence of paths to the request files
    """
  __slots__=('requestfiles','children')
  ##TODO: initialization must load the children into an attribute (probably `children`) listed in _child_seq_attrs
  ##explain in documentation about the `children` attribute
  ##all the dependencies, etc. (how doit knows what is up-to-date)

#Handle command-line execution
if __name__ == '__main__':
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('requestfile', help="Path to file containing the request(s) to run")
  #TODO: allow selecting a subset of the requests?
  cmdline=parser.parse_args()
  requestfile=filepath.Path(cmdline.requestfile,isFile=True)
  
  #Confirm that specified request file exists
  assert requestfile.exists(), "Could not find specified request file %s"%requestfile.fullpath
  
  #Get an interator for the individual requests
  allreqs = TotipotentRequest.all_from_yaml(requestfile.fullpath)
  
  #Process each request
  for req in allreqs:
    req.run()
  
