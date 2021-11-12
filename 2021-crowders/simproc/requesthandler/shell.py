"""Support for requests that execute shell commands."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call
from shutil import copy2

#This package
from . import request
from . import yaml_manager
from . import logging

logger=logging.getLogger(__name__)

class ShellCommandRequestBase(request.Request):
  """Base class for requests that execute shell commands
  
  Subclasses must do the following:
  
  - define _outputfile_attrs and/or _more_outputfiles, either as class or instance attributes
  - define a property attribute cmd_str that provides the shell command to be executed, as a string"""

  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Run the shell command
    call(self.cmd_str,shell=True)

_GeneralShellCommandRequest_props_schema_yaml="""#GeneralShellCommandRequest
name: {type: string}
outfile: {type: pathlike}
errfile:
  anyOf:
    - {type: pathlike}
    - {type: 'null'}
command: {type: string}"""

class GeneralShellCommandRequest(ShellCommandRequestBase):
  """Request for simple shell commands
  
  User-defined attributes:
  
    - command: string representing command to execute
    - outfile: Path to output file
    - errfile: Path to error output file, or None to redirect to `outfile`"""
  _self_task=True
  _config_attrs=['outfile','errfile','command']
  _validation_schema=request.Request.update_schema(_GeneralShellCommandRequest_props_schema_yaml)
  _validation_schema.required=['name','outfile','command']
  _outputfile_attrs=['outfile']
  @property
  def cmd_str(self):
    cmd="%s >'%s' "%(str(self.command),self.renderstr(self.outfile))
    if getattr(self,'errfile',None) is None:
      cmd += "2>&1"
    else:
      cmd += "2>'%s'"%self.errfile
    return cmd

_CommonShellCommandRequest_props_schema_yaml="""#CommonShellCommandRequest
name: {type: string}
commandargs: {type: array}
input_args:
  type: array
  items: {type: integer}
output_args:
  type: array
  items: {type: integer}
"""

class CommonShellCommandRequest(ShellCommandRequestBase):
  """Request for shell commands of the form <command name> <argument list>

  User-defined attributes:

    - commandargs: list of strings containing the command and all arguments
    - input_args: list of indices in the argument list that contain input files
    - output_args: list of indices in the argument list that contain output files"""
  _self_task=True
  _config_attrs=['commandargs','input_args','output_args']
  _validation_schema=request.Request.update_schema(_CommonShellCommandRequest_props_schema_yaml)
  _validation_schema.required=['name','commandargs']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(CommonShellCommandRequest, self).__init__(**kwargs)
    #Store list of input and output files
    self._more_inputfiles=[self.commandargs[idx] for idx in getattr(self,'input_args',[])]
    self._more_outputfiles=[self.commandargs[idx] for idx in getattr(self,'output_args',[])]
  @property
  def cmd_str(self):
    arglist = [self.renderstr(itm) for itm in self.commandargs]
    return ' '.join(arglist)

_CopyFileRequest_props_schema_yaml="""#CopyFileRequest
name: {type: string}
source: {type: pathlike}
destination: {type: pathlike}
"""

class CopyFileRequest(request.Request):
  """Request to copy a single file from one location to another

  User-defined attributes:

    - source = path to the source file
    - destination = path to the destination file"""
  _self_task=True
  _config_attrs=['source','destination']
  _inputfile_attrs=['source']
  _outputfile_attrs=['destination']
  _validation_schema=request.Request.update_schema(_CopyFileRequest_props_schema_yaml)
  _validation_schema.required=['source','destination']
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Do the copy operation
    copy2(self.renderstr(self.source),self.renderstr(self.destination))

#Register for loading from yaml
yaml_manager.register_classes([GeneralShellCommandRequest, CommonShellCommandRequest, CopyFileRequest])
