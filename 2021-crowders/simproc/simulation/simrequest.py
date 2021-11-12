"""Base functionality for simulation requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import sys

#Site packages
import numpy as np
import fenics as fem
from scipy.sparse import csr_matrix

#This package
from ..requesthandler.commandseq import WithCommandsRequest
from ..requesthandler import yaml_manager
from ..requesthandler import locators
from ..requesthandler import nested
from .meshinfo import MeshInfo
from ..postproc.plotseries import PlotSeries
from . import expressions
from ..requesthandler import logging

logger=logging.getLogger(__name__)

#Locators
locators.folder_structure.update(OwnSolutionFile=['solutions',0,1])
locators.folder_structure.update(OtherSolutionFile=['solutions'])

EmptyConditions_schema_yaml="""#EmptyConditions
type: object
properties: {}
required: []
additionalProperties: False
"""
EmptyConditions_schema=yaml_manager.readstring(EmptyConditions_schema_yaml)
EmptyConditions=nested.WithNested(**EmptyConditions_schema)
EmptyConditions.properties=nested.WithNested(**EmptyConditions.properties)

GenericConditions_props_schema_yaml="""#GenericConditions
elementorder: {type: integer}
family: {type: string}
dirichlet: {type: object}
neumann: {type: object}
"""
GenericConditions_props_schema=yaml_manager.readstring(GenericConditions_props_schema_yaml)
GenericConditions=EmptyConditions.get_copy()
GenericConditions.properties=nested.WithNested(**GenericConditions_props_schema)
GenericConditions.required=['elementorder']

_SimulationRequest_props_schema_yaml="""#SimulationRequest
mesh:
  anyOf:
    - {type: 'null'}
    - {type: pathlike}
meshmeta:
  anyOf:
    - {type: 'null'}
    - {type: pathlike}
hasmeshfuncs:
  type: boolean
conditions:
  type: object
dataextraction:
  type: array
loaddata:
  type: array
metadata:
  type: object
solver_parameters:
  type: object
skipsolve: {type: boolean}
meshinfo: {}
"""

class SimulationRequest(WithCommandsRequest):
  """Base class for FEniCS simulations
  
  The `run` method of this class includes only basic functionality:

    - Do the standard pre-run checks
    - Load the mesh
    - Call the method `run_sim`, which is to be defined by derived classes or customization
    - Process the output commands.
  
  Of course, even this simple behavior can be overidden through customization if desired.
  
  User-defined attributes:
  
    - mesh: path to mesh hdf5 file
    - meshmeta: optional, path to mesh metadata yaml file
    - hasmeshfuncs: optional, False to avoid loading mesh functions
    - conditions: dictionary specifying model conditions such as element order, boundary conditions, etc.
        Many simulation request types will process this dictionary into a ``nested.WithNested`` instance called ``conditions_processed``.
    - dataextraction = a sequence of data extraction commands to execute after solving the model

      Each command is a pair (cmdname, arguments), where:

        - cmdname = name of the simulator object's method to call, as a string
        - arguments = dictionary of arguments to the method: {argname: value,...}

    - loaddata = a sequence of data loading commands to execute at simulator initialization
    
      Each command is a triple (attrpath, filepath, fieldtag), where:
      
        - attrpath = the attribute path of the simulator attribute, as a string, to store the loaded data
        - filepath = path to the data file to load, as string, relative to folderstructure.datafolder
        - fieldtag = string identifying the HDF5 field name to load
    
    - metadata = dictionary of metadata about the model, for use in post-processing"""
  _self_task=True
  _inputfile_attrs=['mesh','meshmeta']
  _config_attrs=['mesh','meshmeta','hasmeshfuncs','conditions','dataextraction','loaddata','metadata']
  _validation_schema=WithCommandsRequest.update_schema(_SimulationRequest_props_schema_yaml)
  _validation_schema.required=['name','mesh','conditions']
  _validation_schema.set_nested('properties.conditions',GenericConditions)

  def __init__(self,**kwargs):
    #Initialization from base class
    super(SimulationRequest, self).__init__(**kwargs)
    #Get input and output files from the command sequences
    self.init_command_sequence('loaddata')
    self.init_command_sequence('dataextraction')

  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Load the mesh, unless provided some other way
    self.load_mesh()
    #Do the simulation
    logger.startTimer("simulation",request_name=getattr(self,"name",None))
    self.run_sim()
    logger.stopTimer("simulation",request_name=getattr(self,"name",None))
    self.sim_timer=logger.timers["simulation"]
    #Generate output, if any
    self.process_command_sequence(attrpath='dataextraction',singlefunc=None,positional=False)
    return

  def load_mesh(self):
    "Load the mesh specified by attributes"
    if self.mesh is None:
      assert getattr(self,'meshinfo',None) is not None, "Must provide either MeshInfo or a mesh file to load."
      ##TODO: dependency checking on the mesh won't work for this, of course
    else:
      #Process mesh-related attributes
      meshmeta=getattr(self,'meshmeta',None)
      if meshmeta is not None:
        meshmeta=self.render(meshmeta)
      hasmeshfuncs=getattr(self,'hasmeshfuncs',True)
      #Load
      self.meshinfo=MeshInfo.load(self.render(self.mesh),meshmeta,hasmeshfuncs)
    return

  def run_sim(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'run_sim' method."%str(type(self)))

  def set_solver_parameters(self):
    """Set the parameters for the FEniCS solver object"""
    for k,v in getattr(self,'solver_parameters',{}).items():
      nested.set_nested(self.solver.parameters, k, v)

  def loadfield(self,attrpath,infpath,fieldtag,idx=None):
    """Load data into the simulator from an HDF5 input file
    
    Arguments:
    
      - attrpath = attribute path to load the data into, as string
      
        Note that this attribute must already exist, and be of the proper type to receive the requested data.
        (Call ``newfunction`` first if needed.)
      
      - infpath = path to input file
      
      - fieldtag = string identifying the HDF5 field name to load
      
      - idx = index number (as integer) specifying location within the given attribute, None (default) if not the attribute itself is not a sequence
    
    No return value."""
    inputval = self.get_nested(attrpath)
    if idx is not None:
      inputval = inputval[idx]
    hdf5=fem.HDF5File(self.meshinfo.mesh.mpi_comm(),str(self.render(infpath)),'r')
    hdf5.read(inputval,fieldtag)
    hdf5.close()

  def newfunction(self,attrpath,spaceattr,funcname):
    """Set up a new fenics function

    Arguments:

      - attrpath = attribute path for the new function

      - spaceattr = attribute path for the fenics function space

      - funcname = name to use for the new function, as a string

    No return value."""
    fspace = self.get_nested(spaceattr)
    newf=fem.Function(fspace,name=funcname)
    self.set_nested(attrpath,newf)

  def setconstant(self,attrpath,constval):
    """Set up a fenics constant

    This is intended for use as a ``loaddata`` command.

    Arguments:
    
      - attrpath = path to save the constant to, as string

      - constval = value of the constant, as number or Stored
  
    No return value."""
    self.set_nested(attrpath,fem.Constant(self.get_stored(constval)))
    return

  def loadexpression(self,attrpath,spaceattr,expression,parameters=None):
    """Set up a fenics expression

    This is intended for use as a ``loaddata`` command.

    Arguments:

    - attrpath = path to save the expression to, as string

    - spaceattr = attribute path to the FunctionSpace to use

    - expression = required string (or Stored), the expression to project
    
      Note that the python ``format`` method is called on the string, using the mesh metadata (NOT simulator metadata) as the keyword arguments.
      This allows the expression to reference variables defining the mesh structure, without using FEniCS parameters.
      
    - parameters = optional list of parameter names (or Stored) to use in the expression, empty for no parameters.
      
      Note that the values of these parameters are taken from the Simulator's ``metadata`` attribute (NOT mesh metadata).
    
    No return value."""
    #Get the functionspace
    fspace=getattr(self,spaceattr)
    #Apply mesh metadata to the expression
    exprstr=self.get_stored(expression).format(**self.meshinfo.metadata)
    #Create the parameters dictionary
    if parameters is None:
      in_params=[]
    else:
      in_params=self.get_stored(parameters)
    params={}
    for k in in_params:
      params[k]=self.metadata[k]
    #Create the expression object
    self.set_nested(attrpath,fem.Expression(exprstr,element=fspace.ufl_element(),**params))
    #Done
    return

  def loadcellmapping(self,attrpath,mapping,fieldtype='scalar'):
    """Set up a fenics UserExpression subclass that varies by cell value

    This is intended for use as a ``loaddata`` command.

    Arguments:

    - attrpath = path to save the expression to, as string

    - fieldtype = optional string specifying type of field:
    
      - 'scalar' (default) for a scalar field
      - 'vector' for a vector field
      - 'matrix' for a rank-2 tensor field

    - mapping = required dictionary, the mapping from cell value to function value

    No return value."""
    #For convenience
    spatial_dims=self.meshinfo.mesh.geometry().dim()
    #Instantiate the appropriate expression subclass
    if fieldtype == 'scalar':
      expr = expressions.VaryingScalarByCell(self.meshinfo.cells,mapping,degree=0)
    elif fieldtype == 'vector':
      expr = expressions.VaryingVectorByCell(self.meshinfo.cells,mapping,spatial_dims,degree=0)
    elif fieldtype == 'matrix':
      expr = expressions.VaryingMatrixByCell(self.meshinfo.cells,mapping,spatial_dims,degree=0)
    else:
      raise Exception("Invalid fieldtype: %s"%fieldtype)
    #Store the result
    self.set_nested(attrpath,expr)
    return

  def process_load_commands(self):
    """Process a list of load data commands

    No return value."""
    self.process_command_sequence(attrpath='loaddata',singlefunc=None,positional=False)
    return

  def writefield_outputfiles(self,outfpath,attrpath='soln',idx=None,outname=None):
    """Compute a list of output files generated by writefield"""
    outlist=[outfpath]
    ofp=self.render(outfpath)
    if ofp.suffix.lower()=='.pvd':
      #There will also be .vtu or .pvtu files
      ##TODO: detect if MPI, and compute pvtu files
      vtufname=ofp.stem + "%06d"%0 + ".vtu"
      vtupath=ofp.folder_path / vtufname
      outlist.append(vtupath)
    return outlist

  def writefield(self,outfpath,attrpath='soln',idx=None,outname=None):
    """Write field to VTK or HDF5 file

    Arguments:

      - outfpath = path to output file

        If the filename ends with ".hdf5" (case insensitive), an HDF5 file is created.
        Otherwise, the filetype is selected by FEniCS.

      - attrpath = attribute path to output, as string, defaults to 'soln'
      
      - idx = index number (as integer) of the solution field to write out, None (default) if not a sequence
      
      - outname = optional output field name within output file, as string, supported only for HDF5 format.
      
        If not provided, output field name defaults to value calculated from attrpath and idx.

    Required attributes:

      - the specified attribute containing the FEniCS field

    No new attributes.

    No return value.

    Output file is written."""
    output = self.get_nested(attrpath)
    if idx is not None:
      output = output[idx]
    if str(self.render(outfpath))[-5:].lower()=='.hdf5':
      #HDF5 format
      if outname is None:
        fieldtag=attrpath
        if idx is not None:
          fieldtag+='_%d'%idx
      else:
        fieldtag = outname
      hdf5=fem.HDF5File(self.meshinfo.mesh.mpi_comm(),self.renderstr(outfpath),'w')
      hdf5.write(output,fieldtag)
      hdf5.close()
    else:
      #Format controlled by FEniCS (including VTK files: .pvd, etc.)
      out_file=fem.File(str(self.render(outfpath)))
      out_file << output
    return

  def splitfield(self,namewhole,namesplit):
    """Call the split() method of a solution field.

    Arguments:

      - namewhole = name of the field to split
      - namesplit = attribute name to store the result in

    Required attributes:

      - The attribute specified by namewhole.

    New attributes:

      - The attribute specified by namesplit.

    No return value.

    No other side-effects."""
    setattr(self,namesplit,getattr(self,namewhole).split())

  def domain_volume(self,attrpath='volume',dxname='dx'):
    """Get the entire problem domain volume from integration
    
    Arguments:
    
      - attrpath = optional, attribute path for storing results, as string
      - dxname = optional, name of attribute with fenics domain volume measure, defaults to "dx"
    
    New attribute added/overwritten.
    No return value.
    No output files."""
    dx=getattr(self,dxname)
    volume=fem.assemble(fem.Constant(1)*dx)
    self.set_nested(attrpath,volume)

  def cell_volume(self,pcell,attrpath):
    """Compute the volume of the specified cell

    Arguments:

      - pcell = physical cell number for the cell to calculate area of
      - attrpath = attribute path for storage of result
    
    Required attributes:

      - meshinfo.mesh = FEniCS Mesh object
      - meshinfo.cells = FEniCS MeshFunction object for cell numbers
    
    No return value."""
    this_dx=fem.Measure('cell',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.cells)
    calcvol=fem.assemble(fem.Constant(1)*this_dx(pcell))
    self.set_nested(attrpath,calcvol)
    return

  def cell_integral(self,attrpath,pcell=None,funcpath=None):
    """Compute the specified integral over the given cell.

    Arguments:

      - attrpath = attribute path for storage of result
      - pcell = optional, physical cell number for the cell to calculate integral over,
          empty or None to use entire model
      - funcpath = optional, path to function to integrate,
          empty or None to use a constant value of 1.0

    Required attributes:

      - meshinfo.mesh = FEniCS Mesh object
      - meshinfo.cells = FEniCS MeshFunction object for cell numbers
    
    No return value."""
    if funcpath is None:
      this_func=fem.Constant(1.0)
    else:
      this_func=self.get_nested(funcpath)
    if pcell is None:
      this_dx=fem.Measure('cell',domain=self.meshinfo.mesh)
    else:
      this_measure=fem.Measure('cell',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.cells)
      this_dx=this_measure(pcell)
    result=fem.assemble(this_func*this_dx)
    self.set_nested(attrpath,result)

  def facet_area(self,pfacet,attrpath,internal=False):
    """Compute the area of the specified facet.

    Arguments:

      - pfacet = physical facet number of facet to calculate area of
      - attrpath = attribute path for storage of result
      - internal = boolean, default False, True for internal boundary, False for external

    Required attributes:

      - meshinfo.mesh = FEniCS Mesh object
      - meshinfo.facets = FEniCS MeshFunction object for facet numbers

    No return value."""
    if internal:
      integral_type='interior_facet'
    else:
      integral_type='exterior_facet'
    this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    calcarea=fem.assemble(fem.Constant(1)*this_ds(pfacet))
    self.set_nested(attrpath,calcarea)
    return

  def facet_integral(self,attrpath,pfacet,funcpath=None,internal=False):
    """Compute the specified integral over the given facet.

    Arguments:

      - attrpath = attribute path for storage of result
      - pfacet = physical facet number of facet to calculate integral over
      - funcpath = optional, path to function to integrate,
          empty or None to use a constant value of 1.0
      - internal = boolean, default False, True for internal boundary, False for external

    Required attributes:

      - meshinfo.mesh = FEniCS Mesh object
      - meshinfo.facets = FEniCS MeshFunction object for facet numbers

    No return value."""
    if funcpath is None:
      this_func=fem.Constant(1.0)
    else:
      this_func=self.get_nested(funcpath)
    if internal:
      integral_type='interior_facet'
    else:
      integral_type='exterior_facet'
    this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    result=fem.assemble(this_func*this_ds(pfacet))
    self.set_nested(attrpath,result)
    return

  def get_pointcoords(self,location):
    """Process a location specifier.

    Arguments:

      - location = specifier of location within the mesh

        This should be a tuple, with length matching the problem dimensions.
        Each entry is either a number or a string.
        Numbers represent physical coordinates within the mesh.
        Strings are replaced with the corresponding entry from the mesh metadata dictionary.
          (which is a required attributed in that case)

    Required attributes:

      - mesh_metadata = only required if needed by location specifiers, dictionary of mesh metadata

    Returns:

      - coords = the converted tuple"""
    coords=tuple()
    for v in location:
      if type(v)==str:
        v=self.meshinfo.metadata[v]
      coords+=(v,)
    return coords

  def line_profile(self,startloc,endloc,num,plotpath,label,attrpath='soln',indep=None,idx=None,metadata=None):
    """Get data to plot a result along the specified line at a single point in time
  
    Arguments:
  
      - startloc = argument to get_pointcoords for start of line
      - endloc = argument to get_pointcoords for end of line
      - num = number of sampled points
      - indep = index of the coordinate parameter to use as the independent variable for the plot (zero-based) (omit to use distance from start point)
      - plotpath = attribute path for storing the generated PlotSeries instance
      - label = series label to assign, as string
      - attrpath = attribute path to data source, defaults to 'soln'
      - indep = identifier for independent variable:
          integer 0-d to use that coordinate of the point, or
          None (default) to use distance from the start point
      - idx = index of the solution field to write out, None (default) if not a sequence
      - metadata = other parameters needed to identify the data series
  
    No return value."""
    #Get the object with the data
    vals=self.get_nested(attrpath)
    if idx is not None:
      vals = vals[idx]
    #Get the points for data extraction
    assert len(startloc)==len(endloc), "Start and end locations have different dimensionality"
    startcoords=self.get_pointcoords(startloc)
    endcoords=self.get_pointcoords(endloc)
    start_ends=[itm for itm in zip(startcoords,endcoords)]
    ranges=[np.linspace(start,end,num) for start,end in start_ends]
    points=[t for t in zip(*ranges)]
    #Function to calculate independent variable for a given point
    if indep is None:
      indep_f = lambda pt: np.sqrt(sum([(startcoords[i]-c)**2 for i,c in enumerate(pt)]))
    else:
      indep_f = lambda pt: pt[indep]
    #Extract data points
    llist=[]
    vlist=[]
    for pt in points:
      try:
        vlist.append(vals(*pt))
        llist.append(indep_f(pt))
      except RuntimeError:
        pass #point is not inside mesh; skip
    #Create PlotSeries
    larr=np.array(llist)
    varr=np.array(vlist)
    series=PlotSeries(xvals=larr,yvals=varr,label=label,metadata=metadata)
    #Store data
    self.set_nested(plotpath,series)
    return

  def calc_norm(self,factor=1,attrpath="soln",outattr="soln_norm",normtype="L2"):
    """Have fenics compute the norm of a function (or fenics vector)

    Arguments:

      - factor = optional, normalization constant, defaults to 1 (useful if set to the Stored domain volume)
      - attrpath = optional, attribute path to the function, defaults to "soln"
      - outattr = optional, attribute path to the output norm, defaults to "soln_norm"
      - normtype = optional, norm type as string, defaults to "L2" """
    divisor=self.get_stored(factor)
    targ=self.get_nested(attrpath)
    res=fem.norm(targ,normtype)/divisor
    self.set_nested(outattr,res)
    return

  def get_solver_matrices(self,attr_A,attr_b):
    """Compute the matrix and vector for the linear algebra problem

    The linear algebra problem must be of the form:

    .. math::
    
      A x = b
    
    Arguments:
    
      - attr_A = attribute path for storing the matrix A, in scipy sparse format
      - attr_b = attribute path for storing the vector b, as a numpy 1D array
      
    Required attributes:
    
      - a = bilinear form
      - L = linear form
      - bcs = boundary conditions"""
    matA, matb = fem.assemble_system(self.a, self.L, self.bcs)
    backmatA=fem.as_backend_type(matA).mat()
    backmatb=fem.as_backend_type(matb).vec()
    assert fem.parameters['linear_algebra_backend'] == 'PETSc', "PETSc backend required."
    A=csr_matrix(A.getValuesCSR()[::-1],shape=backmatA.size)
    b=backmatb.getArray()
    self.set_nested(attr_A,A)
    self.set_nested(attr_b,b)
    return

  # def compute_determinant(self,attr_A,attr_det):
  #   """Compute the matrix determinant
    
  #   Arguments:
    
  #     - attr_A = attribute path to the matrix
  #     - attr_det = attribute path for storing the determinant"""
  #   a=self.get_nested(attr_A)
  #   logger.startTimer("determinant")
  #   det = np.linalg.det(A)
  #   logger.stopTimer("determinant")
  #   self.set_nested(attr_det,det)
  #   return

  def subdomain_to_meshfunction(self,meshfunc_value,attr_class,attr_meshfunc="meshinfo.facets",init_kwargs=None):
    """Apply the specified SubDomain-derived class to a MeshFunction

    Arguments:

      - meshfunc_value = value to use for the MeshFunction on this SubDomain
      - attr_class = attribute path to the subclass of SubDomain
      - attr_meshfunc = optional, attribute path to the MeshFunction
          (defaults to ``meshinfo.facets``, for specifying boundary surfaces)
      - init_kwargs = optional, keyword arguments to be passed to the constructor

    """
    #Process default arguments
    if init_kwargs is None:
      init_kwargs={}
    #Get the SubDomain subclass
    the_class=self.get_nested(attr_class)
    #Get the MeshFunction
    the_meshfunc=self.get_nested(attr_meshfunc)
    #Instantiate the SubDomain subclass
    the_subdomain=the_class(**init_kwargs)
    #Mark the SubDomain
    the_subdomain.mark(the_meshfunc,meshfunc_value)
    return

#Register for loading from yaml
yaml_manager.register_classes([SimulationRequest])

