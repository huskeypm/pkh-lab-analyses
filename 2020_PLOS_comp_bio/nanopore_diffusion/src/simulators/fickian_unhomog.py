"""Solve the unhomogenized Fickian diffusion problem
and extract the data needed for post-processing efforts"""

#Standard library
import os
import os.path as osp

#Site packages
import fenics as fem

#Local
import simulator_general

class UnhomogFickianConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with UnhomogFickianSimulator

  Attributes:

    - dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    - D_bulk = bulk diffusion constant"""
  __slots__=['dirichlet','D_bulk']

class UnhomogFickianSimulator(simulator_general.GenericSimulator):
  """Simulator for Unhomogenized Fickian Diffusion

  Additional attributes not inherited from GenericSimulator:

    - conditions = instance of UnhomogFickianConditions
    - V = FEniCS FunctionSpace on the mesh
    - V_vec = FEniCS VectorFunctionSpace on the mesh
    - bcs = FEniCS BCParameters
    - ds = FEniCS Measure for facet boundary conditions
    - c = FEniCS TrialFunction on V
    - v = FEniCS TestFunction on V
    - a = bilinear form in variational problem
    - L = linear form in variational problem"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(UnhomogFickianSimulator, self).__init__(modelparams)

    #Get conditions
    self.conditions=UnhomogFickianConditions(**modelparams.conditions)

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.meshinfo.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.meshinfo.mesh, "CG", self.conditions.elementorder)

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in self.conditions.dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets) ##TODO: specify which facets are Neumann?
    if hasattr(self.conditions,'neumann'):
      raise NotImplementedError

    #Define variational problem
    self.c=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    self.a=fem.dot(fem.grad(self.c),fem.grad(self.v))*fem.dx
    self.L=fem.Constant(0)*self.v*self.ds

  def run(self):
    "Run this simulation."
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bcs)
    return

simulatorclasses={'fickian_unhomog':UnhomogFickianSimulator}
