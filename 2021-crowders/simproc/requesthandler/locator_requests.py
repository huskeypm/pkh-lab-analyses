"""Requests dealing specifically with locators.
These things can't go in the locators module itself without creating a circular dependency."""

#This package
from . import locators
from . import request
from . import yaml_manager
from . import logging

logger=logging.getLogger(__name__)

_DumpFolderStructure_props_schema_yaml="""#DumpFolderStructure
outfile: {type: pathlike}"""

class DumpFolderStructure(request.Request):
  """Write the current folder structure to a yaml output file

  Note: this is executed at initialization.
  Running the request has no further effect.
  This is done in order to assist in debugging the loading of subsequent requests.

  User-Provided Attributes:

    - outfile: path to the file containing the requests, as instance of filepath.Path or Locator"""
  _validation_schema=request.Request.update_schema(_DumpFolderStructure_props_schema_yaml)
  _validation_schema.required=['outfile']
  _outputfile_attrs=['outfile']
  _self_task=True
  def __init__(self,**kwargs):
    #Initialization from base class
    super(DumpFolderStructure, self).__init__(**kwargs)
    #Run at initialization time
    self.execute()
  def execute(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Get the folder structure to output
    d=dict(locators.folder_structure.items())
    #Write to file, without flow style
    yaml_manager.writefile_flow(d,self.render(self.outfile))
  def run(self):
    pass
  # def run(self):
  #   def self.execute()

#Register for loading from yaml
yaml_manager.register_classes([DumpFolderStructure])
