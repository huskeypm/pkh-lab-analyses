"""Run gmsh"""

#This package
from ..requesthandler.request import Request
from ..requesthandler.shell import ShellCommandRequestBase
from ..requesthandler.yaml_manager import register_classes, readstring as readyaml
from ..requesthandler import locators

#Locators
locators.folder_structure.update(geotemplate=['mesh','templates'])
locators.folder_structure.update(geofile=['mesh','output',0,'geo'])
locators.folder_structure.update(mshfile=['mesh','output',0,'msh'])
locators.folder_structure.update(gmsh_outfile=['mesh','output',0,'gmsh_out'])
locators.folder_structure.update(meshmetafile=['mesh','output',0,'metadata'])

_GmshRequest_props_schema_yaml="""#GmshRequest
integer_arg: {type: integer}
geofile: {type: pathlike}
mshfile: {type: pathlike}
gmsh_outfile: {type: pathlike}
meshmetafile: {type: pathlike}"""

class GmshRequest(ShellCommandRequestBase):
  """Run gmsh
  
  User-defined attributes:
  
    - geofile: Path to input .geo file
    - mshfile: Path to output .msh file
    - gmsh_outfile: Path to text file to store gmsh message output
    - meshmetafile: optional, Path to yaml file to store mesh metadata
    - integer_arg: optional, integer to pass to gmsh on the command line, to specify meshing dimension
        defaults to 0, which indicates that the .geo file contains the appropriate `Mesh` command."""
  _self_task=True
  _outputfile_attrs=['mshfile','gmsh_outfile','meshmetafile']
  _inputfile_attrs=['geofile']
  _config_attrs=['geofile','mshfile','gmsh_outfile','meshmetafile','integer_arg']
  _validation_schema=Request.update_schema(_GmshRequest_props_schema_yaml)
  _required_attrs=['name','geofile','mshfile','gmsh_outfile']
  @property
  def cmd_str(self):
    #Integer argument
    int_arg = getattr(self,'integer_arg',None)
    int_arg = 0 if int_arg is None else int_arg
    cmd = "gmsh -%d -format msh2"%int_arg
    #meshmetafile, if provided
    if getattr(self,'meshmetafile',None) is not None:
      cmd += " -setstring meshmetafile '%s'"%self.renderstr(self.meshmetafile)
    #All the other files
    cmd += " -o '%s' '%s' >'%s'"%(self.renderstr(self.mshfile),self.renderstr(self.geofile),self.renderstr(self.gmsh_outfile))
    return cmd

#Register for loading from yaml
register_classes([GmshRequest])
