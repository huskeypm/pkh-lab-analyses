"""Run gmsh, dolfin-convert, and the HDF5 conversion with a single request"""

#This package
from ..requesthandler import locators
from ..requesthandler.yaml_manager import register_classes, readstring as readyaml
from ..requesthandler.request import Request
from . import gmsh_runner
from . import dconv_runner
from . import hdf5_conv

file_extensions_yaml="""#File Extensions
geofile: .geo
mshfile: .msh
gmsh_outfile: .txt
meshmetafile: .yaml
dconv_outfile: .txt
mesh_xml: .xml
facet_xml: _facet_region.xml
cell_xml: _physical_region.xml
mesh_hdf5file: .hdf5
"""
file_extensions=readyaml(file_extensions_yaml)

#Table of information for child requests:
# child request attribute name, request class, request name suffix, list of attributes for kwargs
child_request_control=[
  ('gmsh_request',gmsh_runner.GmshRequest,'.gmsh',
    ['geofile','mshfile','gmsh_outfile','meshmetafile','integer_args']),
  ('dconv_request',dconv_runner.DolfinConvertRequest,'.dconv',
    ['mshfile','mesh_xml','dconv_outfile']),
  ('hdf5_request',hdf5_conv.HDF5ConvertRequest,'.hdf5',
    ['mesh_xml','facet_xml','cell_xml','mesh_hdf5file'])
]

_GeoToHDF5Request_props_schema_yaml="""#GeoToHDF5Request
mesh_stem: {type: string}
integer_arg: {type: integer}
geofile: {type: pathlike}
mshfile: {type: pathlike}
gmsh_outfile: {type: pathlike}
meshmetafile: {type: pathlike}
dconv_outfile: {type: pathlike}
mesh_xml: {type: pathlike}
facet_xml: {type: pathlike}
cell_xml: {type: pathlike}
mesh_hdf5file: {type: pathlike}
gmsh_request: {}
dconv_request: {}
hdf5_request: {}"""

class GeoToHDF5Request(Request):
  """Run gmsh, dolfin-convert, and the HDF5 conversion with a single request
  
  User-defined attributes:

    - name: name of request, as string, required
    - mesh_stem: stem name of mesh, as string, required
    - integer_arg: optional, integer to pass to gmsh on the command line, to specify meshing dimension
        defaults to 0, which indicates that the .geo file contains the appropriate `Mesh` command.

  Optional attributes for file paths (calculated from locators if not provided):

    - geofile: Path to input .geo file
    - mshfile: Path to output .msh file
    - gmsh_outfile: Path to text file to store gmsh message output
    - meshmetafile: optional, Path to yaml file to store mesh metadata
    - dconv_outfile: Path to text file to store dolfin-convert message output
    - mesh_xml = Path to .xml file for the mesh itself
    - facet_xml = Path to .xml file containing facet meshfunction data
    - cell_xml = Path to .xml file containing cell meshfunction data
    - mesh_hdf5file = Path to output .hdf5 file

  Calculated Attributes:

    - gmsh_request: Request that will run gmsh
    - dconv_request: Request that will run dolfin-convert
    - hdf5_request: Request that will convert to hdf5"""
  _self_task = False
  _outputfile_attrs=['mshfile','gmsh_outfile','meshmetafile','dconv_outfile','mesh_xml','facet_xml','cell_xml','mesh_hdf5file']
  _inputfile_attrs=['geofile']
  _allfile_attrs=_inputfile_attrs+_outputfile_attrs
  _child_attrs=['gmsh_request','dconv_request','hdf5_request']
  _validation_schema=Request.update_schema(_GeoToHDF5Request_props_schema_yaml)
  _validation_schema.required=['name','mesh_stem']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(GeoToHDF5Request, self).__init__(**kwargs)
    #Set up locators
    for attrname in self._allfile_attrs:
      if not hasattr(self,attrname):
        ext=file_extensions[attrname]
        the_loc=getattr(locators,attrname)(self.mesh_stem+ext)
        setattr(self,attrname,the_loc)
    #Set up child requests
    for child_attr, rclass, suffix, keylist in child_request_control:
      kwargs=dict([(k,getattr(self,k)) for k in keylist if hasattr(self,k)])
      kwargs['name']=self.name+suffix
      chreq=rclass(**kwargs)
      setattr(self,child_attr,chreq)

#Register for loading from yaml
register_classes([GeoToHDF5Request])

