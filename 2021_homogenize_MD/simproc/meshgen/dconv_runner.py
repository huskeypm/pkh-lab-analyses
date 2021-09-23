"""Run dolfin-convert"""

#This package
from ..requesthandler.request import Request
from ..requesthandler.shell import ShellCommandRequestBase
from ..requesthandler.yaml_manager import register_classes, readstring as readyaml
from ..requesthandler.filepath import Path
from ..requesthandler import locators

#Locators
locators.folder_structure.update(mesh_xml=['mesh','output',0,'xml'])
locators.folder_structure.update(facet_xml=['mesh','output',0,'xml'])
locators.folder_structure.update(cell_xml=['mesh','output',0,'xml'])
locators.folder_structure.update(dconv_outfile=['mesh','output',0,'dconv_out'])

#Constants
DCONV_XML_SUFFIX={'facet_xml':'_facet_region', 'cell_xml':'_physical_region'}

def get_paths_facet_cell(self):
  """Obtain the paths for the facet and cell xml files from the mesh xml file path
  
  Required attributes:
  
    - mesh_xml = Path or locator to the mesh xml file
  
  New attributes:
  
  - facet_xml = Path to .xml file containing facet meshfunction data
  - cell_xml = Path to .xml file containing cell meshfunction data"""
  global DCONV_XML_SUFFIX
  #Compute paths for facet_xml and cell_xml if not provided
  mesh_xml=self.render(self.mesh_xml)
  for attrname in DCONV_XML_SUFFIX.keys():
    if not hasattr(self,attrname):
      filename=mesh_xml.stem+DCONV_XML_SUFFIX[attrname]+mesh_xml.suffix
      setattr(self,attrname,Path(mesh_xml.folder,filename))

_DolfinConvertRequest_props_schema_yaml="""#DolfinConvertRequest
mshfile: {type: pathlike}
mesh_xml: {type: pathlike}
facet_xml: {type: pathlike}
cell_xml: {type: pathlike}
dconv_outfile: {type: pathlike}"""

class DolfinConvertRequest(ShellCommandRequestBase):
  """Run dolfin-convert
  
  User-defined attributes:
  
    - mshfile: Path to input .msh file
    - mesh_xml = Path to output .xml file
       Two other files are also created by dolfin-convert, in the same directory,
       which contain mesh function data.
       This is the path to output file containing the mesh itself.
    - dconv_outfile: Path to text file to store dolfin-convert message output
  
  Calculated attributes:
  
    - facet_xml = Path to .xml file containing facet meshfunction data
    - cell_xml = Path to .xml file containing cell meshfunction data
  """
  _self_task=True
  _inputfile_attrs=['mshfile']
  _outputfile_attrs=['mesh_xml','dconv_outfile','facet_xml','cell_xml']
  _config_attrs=['mshfile','mesh_xml','dconv_outfile']
  _validation_schema=Request.update_schema(_DolfinConvertRequest_props_schema_yaml)
  _validation_schema.required=['name','mshfile','mesh_xml','dconv_outfile']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(DolfinConvertRequest, self).__init__(**kwargs)
    #Compute paths for facet and cell xml file
    get_paths_facet_cell(self)
  @property
  def cmd_str(self):
    return "dolfin-convert %s %s > %s"%(self.renderstr(self.mshfile),self.renderstr(self.mesh_xml),self.renderstr(self.dconv_outfile))

#Register for loading from yaml
register_classes([DolfinConvertRequest])
