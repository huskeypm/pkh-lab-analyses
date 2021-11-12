"""FEniCS mesh support, and DOF support"""

#Standard library
import os.path as osp

#Site packages
import numpy as np
import fenics as fem

#This package
from ..requesthandler import yaml_manager

class MeshInfo:
  """Bunch of mesh-related data

  Attributes:

    - mesh = FEniCS Mesh
    - facets = FEniCS MeshFunction of gmsh Physical Surface number (3D) or Physical Line number (2D)
    - cells = FEniCS MeshFunction of gmsh Physical Volume number (3D) or Physical Surface number (2D)
    - metadata = dictionary of metadata about the mesh, such as parametric locations

  If a mesh is not provided, its dimensions won't be known,
  so the resulting MeshFunctions might not be very useful.

  A note on the terminology used in FEniCS and gmsh:

  |  The FEniCS information below is from page 185-186 of the FEniCS book.
  |  d = number of dimensions in entity,
  |  D = number of dimensions in problem (maximum entity dimension)
  |  D-d = "codimension" of entity
  |  Terms:
  |    D=2, d=1: fenics facet (facet_region xml) = fenics edge = gmsh physical line
  |    D=2, d=2: fenics cell (physical_region xml) = fenics face = gmsh physical surface
  |    D=3, d=2: fenics facet (facet_region xml) = fenics face = gmsh physical surface
  |    D=3, d=3: fenics cell (physical_region xml) = fenics ____ = gmsh physical volume
  |    also, d=0 is a fenics vertex"""

  def __init__(self,mesh=None,facets=None,cells=None,metadata=None):
    if mesh is None:
      self.mesh=fem.Mesh()
    else:
      self.mesh=mesh
    if facets is None:
      self.facets=fem.MeshFunction("size_t", self.mesh, max(0,self.mesh.geometry().dim()-1)) #Facets are of dimension d-1
    else:
      self.facets=facets
    if cells is None:
      self.cells=fem.MeshFunction("size_t", self.mesh, self.mesh.geometry().dim()) #Cells are of dimension d
    else:
      self.cells=cells
    if metadata is None:
      self.metadata={}
    else:
      self.metadata=metadata

  @classmethod
  def load(cls,mesh_hdf5,meshmetafile,loadfuncs=True):
    """Load Mesh and MeshFunctions from HDF5 file, and mesh metadata from yaml
    
    Arguments:
    
      - mesh_hdf5 = path to the hdf5 file for the mesh
      - meshmetafile = path to the mesh metadata yaml file
      - loadfuncs = optional, True (default) to load mesh functions, False otherwise"""
    #Load mesh metadata file, if it exists
    if meshmetafile is None:
      metadata=None
    else:
      metadata=yaml_manager.readfile(str(meshmetafile))
    #Initialize empty mesh
    mesh=fem.Mesh()
    #Open HDF5 file
    hdf5=fem.HDF5File(mesh.mpi_comm(),str(mesh_hdf5),'r')
    #Get the mesh
    hdf5.read(mesh,'mesh',False)
    #Initialize the object
    self=cls(mesh=mesh,metadata=metadata)
    #Load meshfunctions if requested
    if loadfuncs is True:
      hdf5.read(self.facets,'facets')
      hdf5.read(self.cells,'cells')
    #Done
    hdf5.close()
    return self

  def spatial_dims(self):
    """Number of spatial dimensions in the mesh (not any of the mesh functions)"""
    return self.mesh.geometry().dim()

  def coordinates(self):
    """Coordinates of mesh points as an array.

    Note that this is not the same as the coordinates of the degrees of freedom,
    which requires defining a function space."""
    return self.mesh.coordinates()


class DOFInfo:
  """Bunch of data related to the degrees of freedom (DOFs)

  Attributes:

    - meshinfo = MeshInfo instance
    - V = the function space"""

  def __init__(self,meshinfo,V):
    self.meshinfo=meshinfo
    self.V=V

  def coordinates(self):
    """Coordinates of DOF points as an array"""
    return self.V.tabulate_dof_coordinates().reshape(self.V.dim(),self.meshinfo.spatial_dims())

  def boundary_points(self,markerval=-99.0):
    """Separate the dof coordinates into boundary and non-boundary points

    Arguments:

      - markerval = optional function value to use for marking boundaries

        This value must be compatible with the function space.
        For example, if the function space is scalar, this value must be
        a scalar constant.  
  
    Returns:
    
      - boundary_pts = array of DOF coordinates on the model boundary
      - nonbound_pts = array of DOF coordinates not on the model boundary"""
    #Convenience variables
    Ndof_pts=self.V.dim()
    dofcoords=self.coordinates()
    #Boundary marker function, from a throw-away Dirichlet boundary condition
    bcf=fem.Function(self.V)
    bc=fem.DirichletBC(self.V, markerval, fem.DomainBoundary())
    bc.apply(bcf.vector())
    #Array from that function
    bcf_array=bcf.vector().get_local()
    #Split dof coordinates into boundary and non-boundary
    boundary_pts=[]
    nonbound_pts=[]
    for i in range(Ndof_pts):
        pt=dofcoords[i]
        if fem.near(bcf_array[i],markerval):
            boundary_pts.append(pt)
        else:
            nonbound_pts.append(pt)
    #Convert lists to arrays
    boundary_pts=np.array(boundary_pts)
    nonbound_pts=np.array(nonbound_pts)
    #Done
    return boundary_pts, nonbound_pts
