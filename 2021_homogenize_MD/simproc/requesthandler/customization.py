"""Functions and classes relevant for implementing customization"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import sys
import types

#Site packages

#This package
from . import filepath
from . import locators
from . import request
from . import yaml_manager

#Locator for module files
locators.folder_structure.update(modulefile=['modules'])

def load_modules(module_name_list):
  """Load the specified list of modules.
  
  Each module must be in a location already reachable, as with an import statement
  
  Arguments:
  
    - module_name_list = list of module names, as strings
  
  Returns:
  
    - loaded_module_list = list of loaded modules, as modules"""
  loaded_module_list=[]
  for modname in module_name_list:
    if modname in sys.modules.keys():
      loaded_module=sys.modules[modname]
    else:
      loaded_module=importlib.import_module(modname)
    loaded_module_list.append(loaded_module)
  return loaded_module_list

def load_module_from_path(fpath):
  """Load the python module at the requested location
  
  The module is NOT added to ``sys.modules``.
  
  Arguments:
  
    - fpath = path to the module file, as a filepath.Path instance
  
  Returns:
  
    - the loaded module"""
  module_name = fpath.stem
  file_path = fpath.fullpath
  #Method for python between 2 and 3.5
  if sys.version_info.major == 2 or sys.version_info.minor < 5:
    import imp
    module = imp.load_source(module_name, file_path)
  #Method for python >= 3.5
  else:
    import importlib.util
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
  return module

_CustomizableRequest_props_schema_yaml="""#CustomizableRequest
modules:
  anyOf:
    - {type: 'null'}
    - type: array
      items: {type: pathlike}
methods:
  anyOf:
    - {type: 'null'}
    - type: array
      items: {type: string}
    - type: object
      additionalProperties: {type: string}
initializations:
  anyOf:
    - {type: 'null'}
    - {type: object}
extra:
  anyOf:
    - {type: 'null'}
    - {type: object}
_custom_methods:
  anyOf:
    - {type: 'null'}
    - {type: array}"""

class CustomizableRequest(request.Request):
  """A request that can monkey-patch itself
  
  Customization refers to adding additional attributes and methods to a request at run-time.
  This is done by loading additional python code files as modules.
  The modules decide what functions they wish to allow to become methods of a request.
  The requests decide which of those functions they want, and what names to call them.
  The request can also provide data to the module to be used to initialize it. 
  
  User-defined attributes:
  
    - modules = a sequence of module file paths to be imported.
      If the module contains the variable ``request_methods``,
      as a list of functions (not function names), then those functions will be
      available for binding as methods of the request,
      by using the ``methods`` attribute.
      As such, those functions should have a first argument of ``self``,
      which will be the request to which the method has been bound.

    - methods = dictionary {method_name: function_name} or list (method_name=function_name) of methods to bind.
      The function_name must be for a function listed in the ``request_methods`` variable
      of exactly one of the modules listed in ``modules``.
      This is intended to allow a module to define multiple options for a given method,
      with the request selecting the one it wants.

    - initializations = a dictionary {module name: {variable: value}}

      Upon loading the module, the values in this dictionary will be passed as keyword arguments
        to the function `initialize_module`, if present, within the module.
      Modules listed here but not in `modules` are silently ignored.
      Will using different initalization values within a single python process actually work as expected?

    - extra = dictionary {additional request parameters: assigned value}
    
  Calculated attributes:
  
    - _more_inputfiles: the list of module files is added here
    - _custom_methods: list of the names of the custom methods added to the instance"""
  _validation_schema=request.Request.update_schema(_CustomizableRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(CustomizableRequest, self).__init__(**kwargs)
    #Read the customization attributes, allowing any to be missing
    modules=getattr(self,'modules',[])
    methods=getattr(self,'methods',[])
    initializations=getattr(self,'initializations',{})
    extra=getattr(self,'extra',{})
    #Assign extra attributes
    for k,v in getattr(self,'extra',{}).items():
      setattr(self,k,v)
    #Module files are file dependencies
    self._more_inputfiles=getattr(self,'_more_inputfiles',[])
    self._more_inputfiles+=modules
    #Load modules
    function_name_to_module={} #dictionary mapping addable functions to their home modules
    modpath_list=[self.render(modloc) for modloc in modules]
    for modpath in modpath_list:
      #Load module
      themod=load_module_from_path(modpath)
      modname = themod.__name__
      #Intialize, if requested
      kwargs = initializations.get(modname,None)
      if (kwargs is not None) and hasattr(themod,'initialize_module'):
        themod.initialize_module(**kwargs)
      #Track names of addable methods in this module
      for function_obj in getattr(themod,'request_methods',[]):
        function_name_to_module[function_obj.__name__]=themod
    #Convert list of methods into dictionary with key=value
    if isinstance(methods,list):
      methods=dict([(k,k) for k in methods])
    #Bind methods
    self._custom_methods=[]
    for method_name, function_name in methods.items():
      assert function_name in function_name_to_module.keys(), "Function `%s` not listed in `request_methods` variable in any of these modules: %s"%(function_name,[str(m) for m in modpath_list])
      themod = function_name_to_module[function_name]
      function_obj = getattr(themod,function_name)
      assert isinstance(function_obj,types.FunctionType), "%s in module %s is not a function"%(function_name,themod.__name__)
      self._custom_methods.append(method_name)
      setattr(self,method_name,types.MethodType(function_obj,self)) #Must type cast to MethodType in order to get implicit first argument `self`
  def to_dict(self,recursive=True):
    d=super(CustomizableRequest,self).to_dict(recursive)
    #Remove attributes added by customization: these aren't used for validation, configuration, or instantiation
    for attr in getattr(self,'_custom_methods',[]):
      o=d.pop(attr)
    for attr in getattr(self,'extra',{}).keys():
      o=d.pop(attr)
    return d

_PythonPathRequest_props_schema_yaml="""#PythonPathRequest
folders:
  type: array
  items: {type: pathlike}
"""

class PythonPathRequest(request.Request):
  """A request to add some folders to the python path
  
  The folder is added when the request is loaded, so the order of requests matters.

  User-defined attributes:
  
    - folders: a list of folder paths
      Each path is assumed to be relative to the DATAFOLDER, unless it is an absolute path.
  """
  _validation_schema=request.Request.update_schema(_PythonPathRequest_props_schema_yaml)
  _validation_schema.required=['folders']
  _self_task=False #This task actually doesn't do anything when run
  def __init__(self,**kwargs):
    #Initialization from base class
    super(PythonPathRequest, self).__init__(**kwargs)
    #Load the desired paths
    for target in self.folders:
      targpath=self.render(target)
      if not targpath.is_absolute():
        targpath=locators.DATAFOLDER / targpath
      targstr=str(targpath)
      sys.path.append(targstr)
  def run(self):
    """Running this request is a no-op. It is active only at load."""
    pass

_ModuleLoadRequest_props_schema_yaml="""#ModuleLoadRequest
modules:
  type: array
  items: {type: string}
"""

class ModuleLoadRequest(request.Request):
  """A request to load some python modules
  
  This could be used to load modules that define classes loadable from a yaml file.
  (Not the same yaml file where this module is loaded, due to how yaml loading proceeds.)
  Note that the module locations must already be on the python path

  User-defined attributes:
  
    - modules: a list of module names"""
  _validation_schema=request.Request.update_schema(_ModuleLoadRequest_props_schema_yaml)
  _validation_schema.required=['modules']
  _self_task=False #This task actually doesn't do anything when run
  def __init__(self,**kwargs):
    #Initialization from base class
    super(ModuleLoadRequest, self).__init__(**kwargs)
    #Load the desired modules
    load_modules(self.modules)
  def run(self):
    """Running this request is a no-op. It is active only at load."""
    pass

#Register for loading from yaml
yaml_manager.register_classes([PythonPathRequest, ModuleLoadRequest])
