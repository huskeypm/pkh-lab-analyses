"""Convert XML mesh and mesh functions to HDF5 format"""

#This package
from ..requesthandler.request import Request
from ..requesthandler.yaml_manager import register_classes, readstring as readyaml
from ..requesthandler import locators
from .dconv_runner import get_paths_facet_cell
from ..requesthandler import logging

logger=logging.getLogger(__name__)

#Site packages
import fenics as fem

#Locators
locators.folder_structure.update(mesh_hdf5file=['mesh','output',0,'hdf5'])

_HDF5ConvertRequest_props_schema_yaml="""#HDF5ConvertRequest
mesh_xml: {type: pathlike}
facet_xml:
  anyOf:
    - {type: pathlike}
    - {type: 'null'}
cell_xml:
  anyOf:
    - {type: pathlike}
    - {type: 'null'}
mesh_hdf5file: {type: pathlike}"""

class HDF5ConvertRequest(Request):
  """Convert FEniCS Mesh and MeshFunctions from XML to HDF5 format
  
  User-defined attributes:
  
    - mesh_xml = Path to input .xml file
        This is the path to output file containing the mesh itself.
    - facet_xml = optional Path to input .xml file containing facet meshfunction data
        Computed from the mesh_xml attribute if not provided, assuming the same directory
        If the facet xml file is not desired, this attribute must be explicitly set to None.
    - cell_xml = optional Path to input .xml file containing cell meshfunction data
        Computed from the mesh_xml attribute if not provided, assuming the same directory
        If the cell xml file is not desired, this attribute must be explicitly set to None.
    - mesh_hdf5file = Path to output .hdf5 file
  """
  _self_task=True
  _outputfile_attrs=['mesh_hdf5file']
  _inputfile_attrs=['mesh_xml','facet_xml','cell_xml']
  _config_attrs=['mesh_xml','facet_xml','cell_xml','mesh_hdf5file']
  _validation_schema=Request.update_schema(_HDF5ConvertRequest_props_schema_yaml)
  _validation_schema.required=['name','mesh_xml','mesh_hdf5file']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(HDF5ConvertRequest, self).__init__(**kwargs)
    #Compute paths for facet and cell xml file
    get_paths_facet_cell(self)
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Read in the xml files (Mesh and MeshFunctions)
    mesh = fem.Mesh(self.renderstr(self.mesh_xml))
    if self.facet_xml is None:
      facets = None
    else:
      facets = fem.MeshFunction("size_t", mesh, self.renderstr(self.facet_xml))
    if self.cell_xml is None:
      cells = None
    else:
      cells =  fem.MeshFunction("size_t", mesh, self.renderstr(self.cell_xml))
    #Output to HDF5
    hdf5=fem.HDF5File(mesh.mpi_comm(),self.renderstr(self.mesh_hdf5file),'w')
    hdf5.write(mesh,'mesh')
    if facets is not None:
      hdf5.write(facets,'facets')
    if cells is not None:
      hdf5.write(cells,'cells')
    hdf5.close()
    #Done
    return

#Register for loading from yaml
register_classes([HDF5ConvertRequest])
