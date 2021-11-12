"For reading requests from yaml files"

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#This package
from . import filepath
from . import yaml_manager
from . import locators
from . import request
from . import logging

logger=logging.getLogger(__name__)

MULTIDOC_DEFAULT=False

#Locators
locators.folder_structure.update(RequestFile=['requests'])
locators.folder_structure.update(RequestTemplate=['requests','templates'])
locators.folder_structure.update(GeneratedRequest=['requests','generated'])

_RequestFileRequest_props_schema_yaml="""#RequestFileRequest
requestfile: {type: pathlike}
multidoc: {type: boolean}
delay_load: {type: boolean}
skip_nonexist: {type: boolean}
_children: {type: array}"""

class RequestFileRequest(request.Request):
  """Request to run all the requests in a given file
  
  User-Provided Attributes:
  
    - requestfile: path to the file containing the requests, as instance of filepath.Path or Locator
    - multidoc: optional boolean, defaults to value of module variable ``MULTIDOC_DEFAULT``, which is initially False
        If True, the file is a multi-document yaml file, with each document providing a request.
        If False, the file should be a single document, which contains a single list of requests.
    - skip_nonexist: optional boolean, defaults to False
        If True, no error is raised if the request file does not exist.
        If False, an error is raised if the request file does not exist.
    - delay_load: optional boolean, NOT RECOMMENDED, defaults to False
        If False, the child requests are loaded at initialization of this request.
        If True, the child requests are not loaded until this request is run.
        CAUTION: this option is not supported for doit.
        The requests defined in the request file will not be run doit if the loading is delayed.
        This also breaks cleanup requests as well.
    
  Calculated Attributes:
  
    - _children: A list storing all child requests"""
  _validation_schema=request.Request.update_schema(_RequestFileRequest_props_schema_yaml)
  _validation_schema.required=['requestfile']
  _child_seq_attrs=['_children']
  _inputfile_attrs=['requestfile']
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileRequest, self).__init__(**kwargs)
    #Default attributes
    self.multidoc=getattr(self,'multidoc',MULTIDOC_DEFAULT)
    self.delay_load=getattr(self,'delay_load',False)
    self.skip_nonexist=getattr(self,'skip_nonexist',False)
    if self.delay_load:
      self._children=[]
    else:
      self.load_children()
  def load_children(self):
    "Load the requests defined in the file"
    #Confirm that the file exists
    reqfpath=self.render(self.requestfile)
    if reqfpath.exists():
      #Load all objects from the yaml file
      allobj = yaml_manager.readfile(reqfpath,self.multidoc)
      #Store child objects that are Request subclasses
      self._children=[ch for ch in allobj if isinstance(ch,request.Request)]
    else:
      if self.skip_nonexist:
        self._children=[]
      else:
        raise AssertionError("Request file does not exist: %s"%str(reqfpath))
  def run(self):
    logger.info("Running requests from file",request_class=type(self).__name__,request_name=getattr(self,"name",None),source_file=self.renderstr(self.requestfile))
    if self.delay_load:
      self.load_children()
    #Run as defined by parent class
    super(RequestFileRequest,self).run()

_RequestFileListRequest_props_schema_yaml="""#RequestFileListRequest
requestfiles:
  type: array
  items: {type: pathlike}
multidoc: {type: boolean}
delay_load: {type: boolean}
_children: {type: array}"""

class RequestFileListRequest(request.Request):
  """Request to run all of the requests in a given list of files
  
  User-Provided Attributes:
  
    - requestfiles: sequence of paths to the request files, each an instance of filepath.Path or Locator
    - multidoc: optional boolean, defaults to True, sets the multidoc option for each child request
  
  Calculated Attributes:
  
    - _children: A list storing all child requests
    """
  _validation_schema=request.Request.update_schema(_RequestFileListRequest_props_schema_yaml)
  _validation_schema.required=['requestfiles']
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileListRequest, self).__init__(**kwargs)
    self._more_inputfiles=self.requestfiles
    #Each listed file is a RequestFileRequest
    self.multidoc=getattr(self,'multidoc',MULTIDOC_DEFAULT)
    self.delay_load=getattr(self,'delay_load',False)
    self._children=[RequestFileRequest(requestfile=ch,multidoc=self.multidoc,delay_load=self.delay_load) for ch in self.requestfiles]

#Register locators and default folder structure
locators.folder_structure.update(RequestFile=['requests',0,1,2,3])

#Register for loading from yaml
yaml_manager.register_classes([RequestFileRequest, RequestFileListRequest])
