"""Support for requests that fill in templates."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call

#Site packages
from jinja2 import Environment, FileSystemLoader

#This package
from . import yaml_manager
from . import commandseq
from . import nested
from . import logging

logger=logging.getLogger(__name__)

_TemplateFileRequest_props_schema_yaml="""#TemplateFileRequest
name: {type: string}
searchpaths:
  type: array
  items: {type: pathlike}
tmplfile: {type: pathlike}
outfile: {type: pathlike}
data: {type: object}
prepcommands: {type: array}
postcommands: {type: array}
"""

class TemplateFileRequest(commandseq.WithCommandsRequest):
  """General and base class for requests to fill in jinja2 template files
  
  You can subclass this by overriding the `get_template_input` method.
  You can also get the same effect using the customization interface
  to load a replacement `get_template_input` method from a module.
  Otherwise, the template input data is assumed to reside in `data`.
  
  User-defined attributes:
  
    - searchpaths: optional list of folders to add to the jinja2 search path for template inheritance
    - tmplfile: path to the input template file, as Path or string
    - outfile: path to the output file, as Path or string
    - data: dictionary of data used to compute the template input values
    - prepcommands = sequence of commands to execute before template generation (e.g. to load additional data)
    - postcommands = sequence of commands to execute after template generation (e.g. to output additional data)"""
  _self_task=True
  _config_attrs=('searchpaths','tmplfile','outfile','data','modules','initializations','extra')
  _inputfile_attrs=['tmplfile']
  _outputfile_attrs=['outfile']
  _validation_schema=commandseq.WithCommandsRequest.update_schema(_TemplateFileRequest_props_schema_yaml)
  _validation_schema.required=['name','tmplfile','outfile','data']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(TemplateFileRequest, self).__init__(**kwargs)
    #Get input and output files from the command sequences
    self.init_command_sequence('prepcommands')
    self.init_command_sequence('postcommands')
  def get_template_input(self):
    return self.data
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Run prepcommands
    self.process_command_sequence(attrpath='prepcommands',singlefunc=None,positional=False)
    #Add the folder containing the requested template file to the search paths
    tmpl_fpath=self.render(self.tmplfile)
    searchpaths=[str(tmpl_fpath.parent)]
    #Template search paths requested
    raw_searchpaths=getattr(self,'searchpaths',[])
    searchpaths+=[self.renderstr(itm) for itm in raw_searchpaths]
    #Load the template
    env=Environment(loader=FileSystemLoader(searchpaths),extensions=['jinja2.ext.do'],trim_blocks=True,keep_trailing_newline=True)
    tmpl=env.get_template(tmpl_fpath.name)
    #tmpl=Template(tdata,trim_blocks=True,keep_trailing_newline=True)
    #Do the calculations for the template values
    input_data=self.get_template_input()
    #Apply the data to the template
    out_data=tmpl.render(**input_data)
    #Write the output file
    with open(self.renderstr(self.outfile),'w') as fh:
      fh.write(out_data)
    #Run postcommands
    self.process_command_sequence(attrpath='postcommands',singlefunc=None,positional=False)

_LocatorRenderingTemplateFileRequest_props_schema_yaml="""#LocatorRenderingTemplateFileRequest
renderlocs:
  type: array
  items: {type: string}
"""

class LocatorRenderingTemplateFileRequest(TemplateFileRequest):
  """Subclass of TemplateFileRequest that renders locators within the provided data before sending to the template.
  
  User-defined attributes (in addition to those of TemplateFileRequest):
  
    - renderlocs: list of attribute paths (relative to ``data``, not the request) containing locators to render
  """
  _validation_schema=TemplateFileRequest.update_schema(_LocatorRenderingTemplateFileRequest_props_schema_yaml)
  _validation_schema.required=TemplateFileRequest._validation_schema.required + ['renderlocs']
  def get_template_input(self):
    #Go through the list of attribute locations
    for addr in self.renderlocs:
      #Get the rendering of the locator at that location
      newval=self.renderstr(nested.get_nested(self.data,addr))
      #Set the data attribute at that location to the rendered string
      nested.set_nested(self.data,addr,newval)
    #Done
    return self.data

#Register for loading from yaml
yaml_manager.register_classes([TemplateFileRequest,LocatorRenderingTemplateFileRequest])
