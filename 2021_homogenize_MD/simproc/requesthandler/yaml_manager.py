"""Interface for dealing with yaml-loadable classes."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
import io

#Site packages
import ruamel.yaml

#This package

#YAML setup

#to avoid the bug identified in
#https://stackoverflow.com/questions/53396845/with-python-ruamel-yaml-lost-anchor-when-loading-in-round-trip-mode
# #Unfortunately, this didn't fix that bug.
# class MyConstructor(ruamel.yaml.constructor.RoundTripConstructor):
#   def construct_yaml_map(self, node):
#     data = ruamel.yaml.comments.CommentedMap()
#     data._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
#     yield data
#     self.construct_mapping(node, data, deep=True)
#     self.set_collection_style(data, node)

# MyConstructor.add_constructor(u'tag:yaml.org,2002:map', MyConstructor.construct_yaml_map)

#Dictionary of all registered classes, by their names
all_registered=odict()

#Placeholder for tracking file being loaded
now_loading=[]

def register_classes(class_list):
  """Register the classes that might be loaded from a yaml file, from a list of classes
  
  Arguments:
  
    - class_list = sequence of classes, as classes"""
  for yclass in class_list:
    yaml.register_class(yclass)
    all_registered[yclass.__name__]=yclass

def newloader(yfile=None,flowstyle=False):
  """For when you need to load more than one yaml file at a time.
  
  Arguments:
  
    - yfile = optional name of file the loader will be used for.
    - flowstyle = optional boolean, to be used for ``default_flow_style``, (use None for ruamel.yaml default)
  
  The returned object is an instance of YAML,
  which is not actually called a loader, but I can't figure out the actual name."""
  if yfile is not None:
    global now_loading
    now_loading.append(yfile)
  yy=ruamel.yaml.YAML(typ="safe", pure=True)
  # yy=ruamel.yaml.YAML(typ="rt", pure=True)
  yy.default_flow_style = flowstyle
  # yy.Constructor = MyConstructor
  for yclass in all_registered.values():
    yy.register_class(yclass)
  return yy

yaml=newloader()
yamlflow=newloader(flowstyle=None)

def filedone():
  """Call to indicate that loading of a file is complete"""
  global now_loading
  now_loading.pop()

def readstring(s):
  """Syntactic sugar for yaml.load()"""
  return yaml.load(s)

def writestring(obj,yy=yaml):
  """like json.dumps but for yaml"""
  with io.StringIO() as strm:
    yy.dump(obj,strm)
    dat=strm.getvalue()
  return dat

def readfile(fpath,multidoc=False):
  """Read from a file, assuming others may be read at the same time."""
  with open(str(fpath),'r') as fp:
    dat=fp.read()
  yaml=newloader(fpath)
  if multidoc:
    allobj=yaml.load_all(dat)
  else:
    allobj=yaml.load(dat)
  filedone()
  return allobj

def writefile(obj,fpath):
  """Write object to yaml file, overwriting"""
  with open(str(fpath),'w') as fp:
    yaml.dump(obj,fp)
  return

def writefile_flow(obj,fpath):
  """Write object to yaml file, overwriting, in flow style"""
  with open(str(fpath),'w') as fp:
    yamlflow.dump(obj,fp)
  return
