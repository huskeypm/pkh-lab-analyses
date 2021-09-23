"""Solve the linearized Poisson-Boltzmann equation"""

#Standard library

#Site packages
import numpy as np
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.nested import WithNested
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder

_LPBConditions_props_schema_yaml="""#LPBConditions
kappa: {type: number}
"""

class LPBSimulator(simrequest.SimulationRequest):
  """Simulator for linearized Poisson-Boltzmann equation
  
  User-defined attributes:
  
    - """
  
  _validation_schema=simrequest.SimulationRequest.update_schema(_LPBConditions_props_schema_yaml,'properties.conditions.properties')
  _validation_schema.get_nested('properties.conditions.required').append('kappa')
  
  def run_sim(self):

    #For convenience
    self.conditions_processed=WithNested(**self.conditions)
    conditions=self.conditions_processed
    conditions.family=conditions.get_nested_default("family","CG")

    #Properties of problem domain
    self.kappa = conditions.kappa
    
    #Function Spaces and Functions
    #Function spaces
    self.V=fem.FunctionSpace(self.meshinfo.mesh,conditions.family,conditions.elementorder)
    #Trial and Test Functions
    self.phi=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V,name="Phi")

    #Load solution if provided
    self.process_load_commands()

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    if hasattr(conditions,'neumann'):
      raise NotImplementedError

    #Measures and facet normals
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)

    #Equation term dictionary
    self.eqnterms=equationbuilder.EquationTermDict()

    #Weak form
    #kappa term
    ufl_term1=(self.kappa**2)*self.phi*self.v*fem.dx
    if self.kappa is not None:
      self.eqnterms.add('bilinear_1',ufl_term1,bilinear=True)
    #Second term
    ufl_term2=fem.dot(fem.grad(self.phi),fem.grad(self.v))*fem.dx
    self.eqnterms.add('bilinear_2',ufl_term2,bilinear=True)
    #Linear forms
    self.eqnterms.add('linear',fem.Constant(0)*self.v*self.ds,bilinear=False)

    #Problem and Solver
    self.a=self.eqnterms.sumterms(bilinear=True)
    self.L=self.eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

    #Get solver parameters
    self.set_solver_parameters()

    #Solve
    if not getattr(self,'skipsolve',False):
      self.solver.solve()

#Register for loading from yaml
yaml_manager.register_classes([LPBSimulator])
