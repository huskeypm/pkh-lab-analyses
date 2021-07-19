"""Solve the unhomogenized Smoluchowski diffusion problem
and extract the data needed for post-processing efforts"""

#Standard library
import argparse
import math
import os
import os.path as osp
from collections import OrderedDict as odict
import numpy as np

#Site packages
import fenics as fem

#Local
import common
import simulator_general

class LPBConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with LPBSimulator

  Attributes:

    - dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    - debye_length = Debye length"""
  __slots__=['dirichlet','debye_length']

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
    self.V = other.V_scalar
    self.ds = other.ds

    #Get conditions
    self.conditions=LPBConditions(**modelparams.conditions)

    #Properties of problem domain
    self.lambda_D = self.conditions.debye_length

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
    if self.lambda_D is None:
      self.a=term2
    else:
      self.a=(1/self.lambda_D**2)*self.phi*self.v*fem.dx + term2
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

class SUConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with SUSimulator

  User-Provided Attributes:

    - dirichlet = dictionary of Dirichlet boundary conditions:
      {physical facet number: solution value, ...}
    - neumann = dictionary of Neumann boundary conditions:
      {physical facet number: [normal derivative values, None for no condition]}
      Note that the normal derivative is NOT the flux; divide the flux by -D to get the normal derivative.
    - reactive = dictionary of reactive boundary conditions:
      {physical facet number: [species_in, species_out]}
      where species_in and species_out are the symbols of the reactant and product species, respectively.
    - species = sequence of dictionaries, each defining a Species object
    - beta = 1/kBT for the temperature under consideration, in units compatible with q times the potential
    - potential = dictionary defining simulator_run.ModelParameters for electric potential
  
  Note also that the attribute bclist (inherited), contains Dirichlet conditions on c, rather than cbar.
    That is, the code will do the Slotboom transformation on the Dirichlet boundary conditions."""
  __slots__=['dirichlet','neumann','species','beta','potential','reactive']

class SUSimulator(simulator_general.GenericSimulator):
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
    super(SUSimulator, self).__init__(modelparams)

    #Get conditions
    self.conditions=SUConditions(**modelparams.conditions)

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
    
    #Elements and Function space(s)
    ele = fem.FiniteElement('CG',self.meshinfo.mesh.ufl_cell(),self.conditions.elementorder)
    mele = fem.MixedElement([ele]*self.Nspecies)
    self.V = fem.FunctionSpace(self.meshinfo.mesh,mele)
    self.V_scalar=fem.FunctionSpace(self.meshinfo.mesh,'CG',self.conditions.elementorder)
    self.V_vec=fem.VectorFunctionSpace(self.meshinfo.mesh,'CG',self.conditions.elementorder)

    #Trial Functions
    self.u = fem.TrialFunction(self.V)
    cbarlist=fem.split(self.u)

    #Test Functions
    vlist=fem.TestFunctions(self.V)

    #Solution function(s)
    self.soln=fem.Function(self.V)

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

    #Slootboom transformation for dirichlet condition
    def transform_value(cval,phival,beta_z):
      return cval*math.exp(beta_z*phival)

    #Dirichlet boundary conditions
    pot_d=self.potsim.conditions.dirichlet
    self.bcs=[]
    dirichlet=getattr(self.conditions,'dirichlet',{})
    for psurf,vals in dirichlet.items():
      for s,value in enumerate(vals):
        if value is not None:
          transval=transform_value(value,pot_d[psurf],self.conditions.beta*self.species[s].z)
          fspace=self.V.sub(s) if self.Nspecies>1 else self.V
          setvalue=transval if self.Nspecies>1 else (transval,)
          self.bcs.append(fem.DirichletBC(fspace,fem.Constant(setvalue),self.meshinfo.facets,psurf))

    #Neumann boundary conditions
    self.nbcs = {}
    neumann=getattr(self.conditions,'neumann',{})
    for psurf,vals in neumann.items():
      for s,value in enumerate(vals):
        if value is not None:
          if type(value)==int or type(value)==float:
            if value != 0: #Neumann conditions of zero can be omitted from the weak form, to the same effect
              self.nbcs[(psurf,s)]=fem.Constant(value)
          elif type(value)==list:
            exprstr, exprargs = value
            self.nbcs[(psurf,s)]=fem.Expression(exprstr,element=ele,**exprargs)
    
    #Reactive boundary terms
    #just convert symbols to species indices
    self.rbcs = {}
    reactive=getattr(self.conditions,'reactive',{})
    for psurf,pair in reactive.items():
      spair=[self.species_indices[symb] for symb in pair]
      self.rbcs[psurf]=spair

    #Calculate Dbar for each species
    self.Dbar_dict={}
    self.Dbar_proj=[]
    for s,spec in enumerate(self.species):
      Dbar=spec.D*fem.exp(-self.conditions.beta*spec.z*self.potsim.soln)
      self.Dbar_dict[s]=Dbar
      self.Dbar_proj.append(fem.project(Dbar,self.V_scalar,solver_type="cg",preconditioner_type="amg")) #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)

    #Weak Form
    allterms=simulator_general.EquationTermDict(simulator_general.EquationTerm)
    #Body terms
    for s,cbar in enumerate(cbarlist):
      if self.species[s].D is not None:
        termname='body_%d'%s
        allterms.add(termname,self.Dbar_dict[s]*fem.dot(fem.grad(cbar),fem.grad(vlist[s]))*self.dx,bilinear=True)
        #Not knowing if any boundary term will be added here, go ahead and apply zero
        termname='boundary_zero_%d'%s
        allterms.add(termname,fem.Constant(0)*vlist[s]*self.dx,bilinear=False)
    #Boundary terms for Neumann conditions
    for tup,expr in self.nbcs.items():
      psurf,s = tup
      termname='boundary_neumann_%d_%d'%(s,psurf)
      bterm = expr*vlist[s]*self.ds(psurf)
      allterms.add(termname,bterm,bilinear=False)
    #Boundary terms for reactive boundaries
    for psurf,pair in self.rbcs.items():
      r,p=pair
      termname='reactive_%d'%psurf
      bterm = self.Dbar_dict[r]*fem.dot(self.n,fem.grad(cbarlist[r]))*(vlist[r]-vlist[p])*self.ds(psurf)
      allterms.add(termname,-bterm,bilinear=True)

    #Problem and Solver
    self.a=allterms.sumterms(bilinear=True)
    self.L=allterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)
    #Set solver parameters to avoid UMFPACK out-of-memory error
    #iterative solver
    ##self.solver.parameters['linear_solver']='cg' #Conjugate Gradient method, an iterative Krylov solver
    ##self.solver.parameters['preconditioner']='amg' #Algebraic MultiGrid preconditioner
    #mumps solver
    #self.solver.parameters['linear_solver']='mumps' #MUMPS, a parallel LU solver
    #gmres with ilu preconditioner
    self.solver.parameters['linear_solver']='gmres'
    self.solver.parameters['preconditioner']='ilu'

  def run(self):
    "Run this simulation"
    #solve
    self.solver.solve()
    #transform back
    self.solnlist=fem.split(self.soln)
    self.clist=[]
    self.cbarlist=[]
    for s,cbar in enumerate(self.solnlist):
      expr=cbar*fem.exp(-self.conditions.beta*self.species[s].z*self.potsim.soln)
      c=fem.project(expr,self.V_scalar,solver_type="cg",preconditioner_type="amg")
      self.clist.append(c)
      cbar_single=fem.project(cbar,self.V_scalar,solver_type="cg",preconditioner_type="amg")
      self.cbarlist.append(cbar_single)
    return

  def fluxfield(self,filename, solnattr='soln', idx=None, fluxattr='flux', D_bulk=None):
    """Flux as vector field (new attribute, and VTK file)"""
    soln=getattr(self,solnattr)
    if idx is not None:
      soln=soln[idx]
    expr=-self.Dbar_proj[idx]*fem.grad(soln)
    fluxres=fem.project(expr,self.V_vec,solver_type="cg",preconditioner_type="amg") #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)
    setattr(self,fluxattr,fluxres)
    vtk_file=fem.File(osp.join(self.outdir,filename))
    vtk_file << fluxres
    return

  def fluxintegral(self,fluxsurf,name,internal=False,fluxsign=None,normalvar=None,fluxattr='flux',idx=None): ##TODO: store also quadrupled value for unit cell?
    """Flux integral over specified facet

    TODO: update arguments information

    Arguments:
      fluxsurf = physical facet number for flux measurement
      name = name for storage in the results dictionary
      internal = boolean, default False, True to use internal boundary, False for external
      fluxsign = '+' or '-' to specify which direction normal to the facet for flux calculation
        Required only if internal==True
      normalvar = optional variable name to write the facet normal components to, as a sequence
    Required attributes:
      flux = flux as vector field
        This requires a previous call to fluxfield
      mesh = FEniCS Mesh object
      facet = FEniCS MeshFunction object for facet numbers
    No new attributes.
    New item(s) added to results dictionary.
    No return value.
    No output files."""
    n=fem.FacetNormal(self.meshinfo.mesh)
    if internal:
      integral_type='interior_facet'
      assert fluxsign=='+' or fluxsign=='-', "Invalid fluxsign: %s"%str(fluxsign)
      this_n=n(fluxsign)
    else:
      integral_type='exterior_facet'
      this_n=n
    if normalvar is not None:
      self.results[normalvar]=['not_yet_computed'] ##TODO: find a way to get coordinates of the facet normal
    this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    flux=getattr(self,fluxattr)
    if idx is not None:
      flux=flux[idx]
    totflux=fem.assemble(fem.dot(flux,this_n)*this_ds(fluxsurf))
    self.results[name]=totflux
    return

  def effective_D(self,name,totflux_name,area_name,startloc,endloc,attrname='soln',idx=None):
    """Calculate effective diffusion constant

    Arguments:

      - name = name for storage in the results dictionary
      - totflux_name = name of previously calculated total flux in results dictionary

          This requires a previous call to fluxintegral.

      - area_name = name of previously claculated area in results dictionary.

          This requires a previous call to facet_area.

      - startloc = argument to get_pointcoords for start of line
      - endloc = argument to get_pointcoords for end of line
      - attrname = optional, name of attribute containing concentration solution, as string
      - idx = index of the solution field to write out, None (default) if not a sequence

    Required attributes:

      - results[toflux_name] = result from previous call to fluxintegral

    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    #Get the object with the data
    vals=getattr(self,attrname)
    if idx is not None:
      vals = vals[idx]
    #Get the points for data extraction
    assert len(startloc)==len(endloc), "Start and end locations have different dimensionality"
    startcoords=self.get_pointcoords(startloc)
    endcoords=self.get_pointcoords(endloc)
    #Calculate distance between the two points
    deltas=[p[1]-p[0] for p in zip(startcoords,endcoords)]
    delta_s=np.sqrt(sum([d**2 for d in deltas]))
    #Calculate the change in concentration between the two points
    delta_c=vals(*endcoords)-vals(*startcoords)
    #Calculate diffusion constant
    Deff=float(self.results[totflux_name]/self.results[area_name]*delta_s/delta_c)
    self.results[name]=Deff
    return

simulatorclasses={'smol_reactive_surface':SUSimulator}
