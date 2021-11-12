"""Solve the homogenized Fickian diffusion problem (corrector problem)
and extract the data needed for post-processing efforts"""

#Standard library

#Site packages
import numpy as np
import fenics as fem
import ufl

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.nested import WithNested
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder
from . import periodic_boundaries

BOUNDTOL=1e-6

_HomogFickianConditions_props_schema_yaml="""#HomogFickianConditions
elementorder: {type: integer}
family: {type: string}
dirichlet:
  anyOf:
    - {type: array}
    - {type: object}
boundaries: {type: array}
"""

class HomogFickianSimulator(simrequest.SimulationRequest):
  """Simulator for Homogenized Fickian Diffusion
  
  Isotropy of the input diffusion constant is assumed."""
  
  _validation_schema=simrequest.SimulationRequest.update_schema(
          _HomogFickianConditions_props_schema_yaml,'properties.conditions.properties',False)
  _validation_schema.set_nested('properties.conditions.required',['elementorder','boundaries'])
  
  def run_sim(self):

    #For convenience
    self.conditions_processed=WithNested(**self.conditions)
    conditions=self.conditions_processed
    conditions.family=conditions.get_nested_default("family","P")
    
    #Mesh calculations
    spatial_dims=self.meshinfo.mesh.geometry().dim()
    meshcoords=self.meshinfo.mesh.coordinates()
    lowerlims=tuple([np.amin(meshcoords[:,i]) for i in range(spatial_dims)])
    upperlims=tuple([np.amax(meshcoords[:,i]) for i in range(spatial_dims)])
    pairedlims=list(zip(lowerlims,upperlims))

    if hasattr(conditions, 'dirichlet'):
      #The dirichlet conditions are a substitute for periodic boundary conditions
      using_periodic=False
    else:
      #Use periodic boundary conditions
      using_periodic=True
      self.bcs=[]
      if spatial_dims==2:
        pbc = periodic_boundaries.PeriodicBoundary2D(*pairedlims)
      elif spatial_dims==3:
        raise NotImplementedError("Sorry, periodic 3D BCs aren't working yet.")
        #pbc = periodic_boundaries.PeriodicBoundary3D(*pairedlims)
      else:
        raise NotImplementedError("You asked for a simulation on a mesh with %d spatial dimensions."%spatial_dims)

    #Function Spaces and Functions
    if using_periodic:
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, conditions.family, conditions.elementorder, constrained_domain=pbc)
    else:
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, conditions.family, conditions.elementorder)
    self.scalar_V = fem.FunctionSpace(self.meshinfo.mesh, conditions.family, conditions.elementorder)
    #Trial and Test Functions
    chi=fem.TrialFunction(self.V)
    v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V,name='chi')

    if hasattr(conditions,'dirichlet'):
      #Use Dirichlet boundary conditions instead of truly periodic boundary conditions
      if isinstance(conditions.dirichlet,dict):
        self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]
      else:
        #This is a workaround for the meshes that don't have meshfunctions
        val=conditions.dirichlet
        def boundary(x, on_boundary):
          ans = False
          for i in range(spatial_dims):
            for j in range(2):
              ans = ans or fem.near(x[i],pairedlims[i][j],BOUNDTOL)
          ans = ans and on_boundary
          return ans
        self.bcs=[fem.DirichletBC(self.V,val,boundary)]

    #Load diffusion constant
    if hasattr(self,'loaddata'):
      self.D=fem.Function(self.scalar_V)
      self.process_load_commands()
    
    #The index objects
    i=ufl.i
    j=ufl.j

    #Measure and normal for external boundaries
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Equation term dictionary
    eqnterms=equationbuilder.EquationTermDict()

    #Bilinear boundary terms
    for psurf in conditions.boundaries:
      termname="bilinear_boundary_%d"%psurf
      form=self.D*self.n[i]*chi[j].dx(i)*v[j]*self.ds(psurf)
      eqnterms.add(termname,form,bilinear=True)

    #Bilinear body terms
    termname="bilinear_body"
    form=-self.D*chi[j].dx(i)*v[j].dx(i)*self.dx
    eqnterms.add(termname,form,bilinear=True)

    #Linear boundary terms
    for psurf in conditions.boundaries:
      termname="linear_boundary_%d"%psurf
      form=self.D*self.n[i]*v[i]*self.ds(psurf)
      eqnterms.add(termname,form,bilinear=False)

    #Linear body terms
    termname="linear_body"
    form=-self.D*v[i].dx(i)*self.dx
    eqnterms.add(termname,form,bilinear=False)

    #FEniCS Problem and Solver
    self.a=eqnterms.sumterms(bilinear=True)
    self.L=eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

    #Get solver parameters
    self.set_solver_parameters()

    #Solve
    if not getattr(self,'skipsolve',False):
      self.solver.solve()

  def macroscale_diffusion(self,respath="D_macro",attrpath="soln",volpath="volume"):
    """Perform the integral to obtain the homogenized diffusion constant
    
    Isotropy of the input D is assumed, but the output D may be anisotropic or even non-diagonal.

    Arguments:

      - respath = optional, attribute path for storing the result
      - attrpath = optional, attribute path storing the solution (the result for chi)
      - volpath = optional, attribute path storing the unit cell volume.
          
    Required attributes (other than those from simulator_general):
    
      - the attribute given by attrpath
      - dx = FEniCS Measure for cells
      - D = FEniCS Function with the diffusion constant
    
    New attribute created/overwitten.
    No return value.
    No output files."""
    d=self.meshinfo.mesh.geometry().dim()
    volume=self.get_nested(volpath)
    kdelta = lambda i,j: 1 if i==j else 0 #Kronecker delta
    soln=self.get_nested(attrpath)
    gradchi=fem.grad(soln)
    matr=[]
    for ii in range(d):
      row=[]
      for jj in range(d):
        term1=kdelta(ii,jj)*fem.assemble(self.D*self.dx)
        term2=fem.assemble(self.D*gradchi[jj,ii]*self.dx)
        val=(term1-term2)/volume
        row.append(float(val))
      matr.append(row)
    self.set_nested(respath,matr)

#Register for loading from yaml
yaml_manager.register_classes([HomogFickianSimulator])
