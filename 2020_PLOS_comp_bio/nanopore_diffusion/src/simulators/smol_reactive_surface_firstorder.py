"""Solve the unhomogenized Smoluchowski diffusion problem,
as a first-order system,
and extract the data needed for post-processing efforts"""

#Standard library
import argparse
import math
import os
import os.path as osp
from collections import OrderedDict as odict

#Site packages
import fenics as fem

#Local
import common
import simulator_general

class LPBConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with LPBSimulator

  Attributes:

    - dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    - kappa = inverse of the Debye length"""
  __slots__=['dirichlet','kappa']

class LPBSimulator(simulator_general.GenericSimulator):
  """Simulator for linearized Poisson-Boltzmann equation.

  Additional attributes not inherited from GenericSimulator:

    - conditions = instance of LPBConditions
    - lambda_D = Debye length
    - V = FEniCS FunctionSpace on the mesh
    - bcs = list of FEniCS DirichletBC instances
    - ds = FEniCS Measure for facet boundary conditions
    - phi = FEniCS TrialFunction on V
    - v = FEniCS TestFunction on V
    - a = bilinear form in variational problem
    - L = linear form in variational problem"""
  def __init__(self,modelparams,other):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance
      - other = simulator to get mesh from"""

    #Load parameters, init output, mesh setup
    super(LPBSimulator, self).__init__(modelparams)

    #Load mesh and meshfunctions
    self.meshinfo=other.meshinfo
    self.V = other.V_one_scalar
    self.ds = other.ds

    #Get conditions
    self.conditions=LPBConditions(**modelparams.conditions)

    #Properties of problem domain
    self.kappa = self.conditions.kappa

    #Function space for scalars and vectors
    ##self.V = fem.FunctionSpace(self.meshinfo.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in self.conditions.dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    ##self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets) ##TODO: specify which facets are Neumann?
    if hasattr(self.conditions,'neumann'):
      raise NotImplementedError

    #Define variational problem
    self.phi=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    term2=fem.dot(fem.grad(self.phi),fem.grad(self.v))*fem.dx
    if self.kappa is None:
      self.a=term2
    else:
      self.a=(self.kappa**2)*self.phi*self.v*fem.dx + term2
    self.L=fem.Constant(0)*self.v*self.ds
    self.soln=fem.Function(self.V)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)
    #Set solver parameters to avoid UMFPACK out-of-memory error
    #iterative solver
    ##self.solver.parameters['linear_solver']='cg' #Conjugate Gradient method, an iterative Krylov solver
    ##self.solver.parameters['preconditioner']='amg' #Algebraic MultiGrid preconditioner
    #mumps
    self.solver.parameters['linear_solver']='mumps' #MUMPS, a parallel LU solver

  def run(self):
    "Run this simulation"
    self.solver.solve()
    return

#Lookup of electric potential simulators by name
potentialsimulatorclasses={'linear_pb':LPBSimulator}

class Species(common.ParameterSet):
  """Information for a single chemical species.
  
  User-Provided attributes:
  
    - symbol: chemical symbol, as string
    - z: ionic charge, as number
    - D: diffusion constant, as number"""
  __slots__=('symbol','z','D')
  _required_attrs=__slots__

class SUFOConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with SUSimulator

  User-Provided Attributes:

    - solver_parameters = dictionary of solver parameters
    - concentrationBC = dictionary of concentration boundary conditions:
      {physical facet number: [concentration values, None for no condition] ...}
    - fluxBC = dictionary of flux boundary conditions:
      {physical facet number: [flux values, None for no condition], ...}
      Note that the normal derivative is NOT the flux; divide the flux by -D to get the normal derivative.
    - reactiveBC = dictionary of reactive boundary conditions:
      {physical facet number: [species_in, species_out]}
      where species_in and species_out are the symbols of the reactant and product species, respectively.
    - species = sequence of dictionaries, each defining a Species object
    - beta = 1/kBT for the temperature under consideration, in units compatible with q times the potential
    - potential = dictionary defining simulator_run.ModelParameters for electric potential
  
  Note also that the concentration boundary conditions are specified as c, rather than cbar.
    That is, the code will do the necessary Slotboom transformation on the concentration boundary conditions.
    The Slotboom transformation does not affect the fluxes."""
  __slots__=['solver_parameters','concentrationBC','fluxBC','reactiveBC','species','beta','potential']

class SUFOSimulator(simulator_general.GenericSimulator):
  """Simulator for Unhomogenized Smoluchowski Diffusion

  Additional attributes not inherited from GenericSimulator:

    - conditions = instance of SUConditions
    - V = FEniCS FunctionSpace on the mesh
    - V_vec = FEniCS VectorFunctionSpace on the mesh
    - bcs = FEniCS BCParameters
    - ds = FEniCS Measure for facet boundary conditions
    - potsim = instance of simulator for the electric potential
    - Dbar = FEniCS Function
    - cbar = FEniCS TrialFunction on V
    - v = FEniCS TestFunction on V
    - a = bilinear form in variational problem
    - L = linear form in variational problem"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(SUFOSimulator, self).__init__(modelparams)

    #Spatial dimensions
    num_dim=self.meshinfo.mesh.geometry().dim()

    #Get conditions
    self.conditions=SUFOConditions(**modelparams.conditions)

    #Species
    self.species=[]
    self.species_dict=odict()
    self.species_indices=odict()
    for s,d in enumerate(self.conditions.species):
      spec=Species(**d)
      self.species.append(spec)
      self.species_dict[spec.symbol]=spec
      self.species_indices[spec.symbol]=s
    self.Nspecies=len(self.species)
    self.Nunk=2*self.Nspecies #2 unknowns per species: flux and concentration
    
    #Mixed Finite Elements
    ele_scalar=fem.FiniteElement('CG',self.meshinfo.mesh.ufl_cell(),self.conditions.elementorder)
    ele_vector=fem.VectorElement('CG',self.meshinfo.mesh.ufl_cell(),self.conditions.elementorder)
    mele=fem.MixedElement([ele_scalar]*self.Nspecies+[ele_vector]*self.Nspecies)

    #Function spaces
    self.V_one_scalar=fem.FunctionSpace(self.meshinfo.mesh,ele_scalar)
    self.V_one_vector=fem.FunctionSpace(self.meshinfo.mesh,ele_vector)
    self.V_all=fem.FunctionSpace(self.meshinfo.mesh,mele)

    #Trial Functions
    u = fem.TrialFunction(self.V_all)
    ulist = fem.split(u)
    ss_cbarlist = ulist[0:self.Nspecies]
    ss_jlist = ulist[self.Nspecies:self.Nunk]

    #Test Functions
    vlist = fem.TestFunctions(self.V_all)
    ss_nulist = vlist[0:self.Nspecies]
    ss_taulist = vlist[self.Nspecies:self.Nunk]

    #Solution function(s)
    self.soln=fem.Function(self.V_all)

    #Measure and normal for external boundaries
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n=fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Set up electric potential field
    potentialparams_dict=self.conditions.potential
    for key in ['modelname','meshname','basename']:
      potentialparams_dict[key]=getattr(modelparams,key)
    potentialparams=simulator_general.ModelParametersBase(**potentialparams_dict)
    self.potsim=potentialsimulatorclasses[potentialparams.equation](potentialparams,self)
    self.potsim.diskwrite=False
    self.potsim.run()
    self.potsim.create_output()
    self.info['potential']=self.potsim.info
    self.outdata.plots=self.potsim.outdata.plots

    #Slootboom transformation for concentration boundary conditions
    def transform_value(cval,phival,beta_z):
      return cval*math.exp(beta_z*phival)

    #Initialze list of essential boundary conditions
    self.bcs=[]

    #cbar boundary conditions
    pot_d=self.potsim.conditions.dirichlet
    concBC=getattr(self.conditions,'concentrationBC',{})
    for psurf,vals in concBC.items():
      for s,value in enumerate(vals):
        if value is not None:
          transval=transform_value(value,pot_d[psurf],self.conditions.beta*self.species[s].z)
          self.bcs.append(fem.DirichletBC(self.V_all.sub(s),fem.Constant(transval),self.meshinfo.facets,psurf))

    #flux boundary conditions
    fluxBC=getattr(self.conditions,'fluxBC',{})
    for psurf,vals in fluxBC.items():
      for s,value in enumerate(vals):
        if value is not None:
          self.bcs.append(fem.DirichletBC(self.V_all.sub(s+self.Nspecies),fem.Constant(value),self.meshinfo.facets,psurf))
    
    #Reactive boundary terms
    #just convert symbols to species indices
    self.rbcs = {}
    reactiveBC=getattr(self.conditions,'reactiveBC',{})
    for psurf,pair in reactiveBC.items():
      spair=[self.species_indices[symb] for symb in pair]
      self.rbcs[psurf]=spair

    #Calculate Dbar for each species
    self.Dbar_dict={}
    self.Dbar_proj=[]
    for s,spec in enumerate(self.species):
      Dbar=spec.D*fem.exp(-self.conditions.beta*spec.z*self.potsim.soln)
      self.Dbar_dict[s]=Dbar
      self.Dbar_proj.append(fem.project(Dbar,self.V_one_scalar,solver_type="cg",preconditioner_type="amg")) #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)

    #Weak Form
    allterms=simulator_general.EquationTermDict(simulator_general.EquationTerm)
    #Body terms
    for s,cbar in enumerate(ss_cbarlist):
      if self.species[s].D is not None:
        #Bilinear Term 1
        termname='body_%d_1'%s
        ufl=ss_nulist[s]*fem.div(ss_jlist[s])*self.dx
        allterms.add(termname,ufl,bilinear=True)
        #Bilinear Term 2
        termname='body_%d_2'%s
        ufl=fem.dot(ss_taulist[s],ss_jlist[s])*self.dx
        allterms.add(termname,ufl,bilinear=True)
        #Bilinear Term 3
        termname='body_%s_3'%s
        ufl=self.Dbar_dict[s]*fem.dot(ss_taulist[s],fem.grad(ss_cbarlist[s]))*self.dx
        allterms.add(termname,ufl,bilinear=True)
        #Zero Term 1
        termname='body_zero_%s_1'%s
        ufl=fem.Constant(0)*ss_nulist[s]*self.dx
        allterms.add(termname,ufl,bilinear=False)
        #Zero Term 2
        termname='body_zero_%s_2'%s
        ufl=fem.dot(fem.Constant((0,)*num_dim),ss_taulist[s])*self.dx
        allterms.add(termname,ufl,bilinear=False)
    #Reactive boundaries
    for psurf,pair in self.rbcs.items():
      r,p=pair
      termname='reactive_%d'%psurf
      ufl=fem.dot(ss_taulist[r]+ss_taulist[p],ss_jlist[r]+ss_jlist[p])*self.ds(psurf)
      allterms.add(termname,ufl,bilinear=True)

    #Problem and Solver
    self.a=allterms.sumterms(bilinear=True)
    self.L=allterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)
    #Get solver parameters from the conditions
    for k,v in getattr(self.conditions,'solver_parameters',{}).items():
      self.solver.parameters[k]=v

  def run(self):
    "Run this simulation"
    #solve
    self.solver.solve()
    #Split, transform, and store
    solnlist=fem.split(self.soln)
    unproj_cbarlist=solnlist[0:self.Nspecies]
    unproj_fluxlist=solnlist[self.Nspecies:self.Nunk]
    self.cbarlist=[]
    self.clist=[]
    self.fluxlist=[]
    for s,cbar in enumerate(unproj_cbarlist):
      #cbar
      cbar_single=fem.project(cbar,self.V_one_scalar,solver_type="richardson",preconditioner_type="jacobi")
      self.cbarlist.append(cbar_single)
      #c
      expr=cbar*fem.exp(-self.conditions.beta*self.species[s].z*self.potsim.soln)
      c=fem.project(expr,self.V_one_scalar,solver_type="cg",preconditioner_type="amg")
      self.clist.append(c)
    for s,flux in enumerate(unproj_fluxlist):
      flux_single=fem.project(flux,self.V_one_vector,solver_type="cg",preconditioner_type="amg")
      self.fluxlist.append(flux_single)

simulatorclasses={'smol_reactive_surface_firstorder':SUFOSimulator}
