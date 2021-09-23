"""Solve the multi-species unhomogenized Fickian diffusion problem,
and extract the data needed for post-processing efforts.
A single species is assumed, as there is no interaction or potential.
Isotropy is assumed, but the diffusion constant may vary spatially."""

#Standard library

#Site packages
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.nested import WithNested
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder
from . import common_methods

class FLSimulator(simrequest.SimulationRequest):
  """Simulator for for Unhomogenized Fickian Diffusion"""

  #Common methods
  calcflux = common_methods.calcflux
  fluxintegral = common_methods.fluxintegral
  effective_D = common_methods.effective_D

  def run_sim(self):

    #For convenience
    self.conditions_processed=WithNested(**self.conditions)
    conditions=self.conditions_processed

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.meshinfo.mesh,'CG',conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.meshinfo.mesh, "CG", conditions.elementorder)

    #Trial Function
    self.c = fem.TrialFunction(self.V)

    #Test Function
    v=fem.TestFunction(self.V)

    #Solution function
    self.soln=fem.Function(self.V,name="soln")

    #Measure and normal for external boundaries
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n=fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Load diffusion coefficient Function
    self.Dlocal=fem.Function(self.V,name="Dlocal")
    self.process_load_commands()

    #Dirichlet boundary conditions
    dirichlet=conditions.get_nested_default("dirichlet",{})
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in dirichlet.items()]

    #Neumann boundary conditions
    self.nbcs = {}
    neumann=getattr(conditions,'neumann',{})
    for psurf,value in neumann.items():
      if value is not None:
        if type(value)==int or type(value)==float:
          if value != 0: #Neumann conditions of zero can be omitted from the weak form, to the same effect
            self.nbcs[psurf]=fem.Constant(value)
        elif type(value)==list:
          exprstr, exprargs = value
          self.nbcs[psurf]=fem.Expression(exprstr,element=ele,**exprargs)
    
    #Weak Form
    allterms=equationbuilder.EquationTermDict()
    #Body term
    termname='body'
    allterms.add(termname,-self.Dlocal*fem.dot(fem.grad(self.c),fem.grad(v))*self.dx,bilinear=True)
    #If no boundary terms will be added, go ahead and apply zero
    if len(self.nbcs)==0:
      termname='boundary_zero'
      allterms.add(termname,fem.Constant(0)*v*self.dx,bilinear=False)
    #Boundary terms for Neumann conditions
    for psurf,expr in self.nbcs.items():
      termname='boundary_neumann_%d'%(psurf)
      bterm = self.Dlocal*expr*v*self.ds(psurf)
      allterms.add(termname,bterm,bilinear=False)

    #Problem and Solver
    self.a=allterms.sumterms(bilinear=True)
    self.L=allterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)
    
    #Get solver parameters
    self.set_solver_parameters()

    #solve
    if not getattr(self,'skipsolve',False):
      self.solver.solve()

#Register for loading from yaml
yaml_manager.register_classes([FLSimulator])




