"""Solve the time-domain unhomogenized PNP diffusion problem
and extract the data needed for post-processing efforts"""

#Standard library

#Site packages
import fenics as fem

#Local
import unitsystem as UN
import common
import simulator_general

class SpeciesInfo(common.ParameterSet):
  """Information on each species

  Attributes:

    - symbol: [list of chemical symbols as strings],
    - z: [list of ionic charges as numbers],
    - initconc: [list of initial concentrations as numbers],
    - D: [list of diffusion constants]
    - N: number of species (calculated)

  Each list must have 1 entry per diffusing chemical species"""
  __slots__=['symbol','z','initconc','D','N']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(SpeciesInfo, self).__init__(**kwargs)
    #Check number of species
    nspec_all=[len(l) for l in kwargs.values()]
    assert min(nspec_all)==max(nspec_all), "Inconsistent number of species: %s"%str(nspec_all)
    self.N=nspec_all[0]

class ReactionInfo(common.ParameterSet):
  """Information on each reaction

  Attributes:

    - constants: [list of reaction rate constants],
    - functions: [list of reaction rate function ....] (all functions must be assigned as methods of the simulator through `customizations`)
    - stoichio: [list of stoichiometric coefficients lists, negative for reactants, positive for products]

      Each entry for stoichiometric coefficients is itself a list (one such list for each reaction),
      with one number for each species in the list for each reaction.
      Thus, the total number of stoichiometric coefficients is the product of the number of reactions and the number of species.
    - N: number of reactions (calculated)

  Each list must have 1 entry per uni-directional reaction (bidirectional reactions are considered as 2 uni-directional reactions each)"""
  __slots__=['constants','functions','stoichio','N']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(ReactionInfo, self).__init__(**kwargs)
    #Check number of reactions
    nreac_all=[len(l) for l in kwargs.values()]
    assert min(nreac_all)==max(nreac_all), "Inconsistent number of reactions: %s"%str(nreac_all)
    self.N=nreac_all[0]

class StoppingCriterion:
  """Define stopping criterion for time-domain simulation.

  Possible attributes:

    - numsteps = number of steps to stop after
    - t_end = time to stop after"""
  def __init__(self,**kwargs):
    assert len(kwargs.keys())==1, "Must specify exactly one stopping criterion, got %d in args: %s"%(len(kwargs.keys()),str(kwargs))
    self.numsteps=kwargs.get('numsteps',None)
    self.t_end=kwargs.get('t_end',None)

class TimeDomainInfo(common.ParameterSet):
  """"Specify time steps and stopping criterion

  Attributes:

    - stepsize: Specification of time step size, as float
    - stopping: instance of StoppingCriterion"""
  __slots__=['stepsize','stopping']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(TimeDomainInfo, self).__init__(**kwargs)
    #Stopping Criterion
    self.stopping=StoppingCriterion(**kwargs['stopping'])

class TDPNPUConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with TDPNPUSimulator

  Attributes:

    - elementorder = see base class
    - dirichlet = dictionary of Dirichlet boundary conditions:
      {physical facet number: [solution values, None for no condition]}
    - neumann = dictionary of Neumann boundary conditions:
      {physical facet number: [normal derivative values, None for no condition]}
      Note that the normal derivative is NOT the flux; divide the flux by -D to get the normal derivative.
    - temperature = the temperature under consideration, as a number
    - eps_r = relative permittivity of the medium
    - species_info = dictionary defining a SpeciesInfo object
    - reaction_info = dictionary defining a ReactionInfo object
    - initial_potential = initial electric potential, assumed constant over space, as number
    - timedomain = instance of simulator_general.TimeDomainInfo
    - beta = optional, calculated from temperature if not provided"""
  __slots__=['dirichlet','neumann','beta','temperature','eps_r','species_info','reaction_info','initial_potential','timedomain']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(TDPNPUConditions, self).__init__(**kwargs)
    #If beta not provided, calculate from temperature
    if not hasattr(self,'beta'):
      self.beta = 1.0/(self.temperature*UN.kB)
    #Get TimeDomainInfo instance
    self.timedomain=TimeDomainInfo(**self.timedomain)

class TDPNPUSimulator(simulator_general.GenericSimulator):
  """Simulator for Unhomogenized Time-Domain Poisson-Nernst-Planck Diffusion

  Additional attributes not inherited from GenericSimulator:

    - conditions = instance of TDPNPUConditions
    - species = instance of SpeciesInfo
    - reactions = instance of ReactionInfo
    - Nspecies = number of chemical species
    - Nvars = number of field variables to solve for
    - dt = timestep
    - V = FEniCS FunctionSpace on the mesh
    - u = FEniCS Function on the FunctionSpace for the current timestep
    - u_k = FENiCS Function on the FunctionSpace for the previous timestep
    - timedeps = time-dependent Expressions that need to be updated at each timestep
    - bcs = FEniCS BCParameters (list of DirichletBC instances)
    - nbcs = dictionary of Neumann boundary conditions: {(facet number, variable index): Expression or Constant for the condition, ...}
    - hbcs = dictionary of PNP hybrid boundary term info : {(): } ##TODO
    - ds = FEniCS Measure for facet boundary conditions
    - n = FEniCS FacetNormal for facet boundary conditions
    - FF = symbolic functional form, which is set equal to zero in the weak form equation
    - J = symbolic Jacobian of FF
    - k = current step number, as integer"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(TDPNPUSimulator, self).__init__(modelparams)

    #Get conditions
    self.conditions=TDPNPUConditions(**modelparams.conditions)
    self.species=SpeciesInfo(**self.conditions.species_info)
    self.reactions=ReactionInfo(**self.conditions.reaction_info)

    #List and count the degrees of freedom
    self.Nspecies=len(self.species.symbol)
    non_species_vars=['Phi']
    varlist=self.species.symbol+non_species_vars
    self.Nvars=len(varlist)

    #Need to keep track of all time-dependent expressions
    self.timedeps=[]

    #Elements and Function space(s)
    ele = fem.FiniteElement('P',self.meshinfo.mesh.ufl_cell(),self.conditions.elementorder)
    mele = fem.MixedElement([ele]*self.Nvars)
    self.V = fem.FunctionSpace(self.meshinfo.mesh,mele)

    #Test and trial functions
    self.u = fem.Function(self.V)
    trialfuncs=fem.split(self.u)
    clist=trialfuncs[:self.Nspecies]
    Phi=trialfuncs[self.Nspecies]
    testfuncs=fem.TestFunctions(self.V)
    vlist=testfuncs[:self.Nspecies]
    vPhi=testfuncs[self.Nspecies]

    #Measure and normal for external boundaries
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n=fem.FacetNormal(self.meshinfo.mesh)

    #Dirichlet boundary conditions
    self.bcs=[]
    dirichlet=getattr(self.conditions,'dirichlet',{})
    for psurf,vals in dirichlet.items():
      for i,value in enumerate(vals):
        if value is not None:
          self.bcs.append(fem.DirichletBC(self.V.sub(i),fem.Constant(value),self.meshinfo.facets,psurf))

    #Neumann boundary conditions
    self.nbcs = {}
    neumann=getattr(self.conditions,'neumann',{})
    for psurf,vals in neumann.items():
      for i,value in enumerate(vals):
        if value is not None:
          if type(value)==int or type(value)==float:
            if value != 0: #Neumann conditions of zero can be omitted from the weak form, to the same effect
              self.nbcs[(psurf,i)]=fem.Constant(value)
          elif type(value)==list:
            exprstr, exprargs = value
            self.nbcs[(psurf,i)]=fem.Expression(exprstr,element=ele,**exprargs)
            if 't' in exprargs.keys():
              self.timedeps.append(self.nbcs[(psurf,i)])

    #Data for hybrid boundary term (hybrid of potential and species)
    self.hbcs = {}
    #First requirement: only surfaces where the potential has either a nonzero Neumann condition, OR any Dirichlet condition
    hbcs_psurfs=[psurf for psurf,vals in neumann.items() if (vals[-1] is not None and vals[-1] != 0)]
    hbcs_psurfs+=[psurf for psurf in dirichlet.keys()]
    #Next requirement: only combinations of surface and species where the species has NO Dirichlet condition
    #BOTH this and the first requirement must be met in order for the term to be nonzero
    for psurf in hbcs_psurfs:
      for s in range(self.Nspecies):
        #Finally, don't include this term for any species that aren't allowed to diffuse
        if self.species.D[s] is not None:
          if psurf not in dirichlet.keys() or dirichlet[psurf][s] is None: #"if there is NOT a Dirichlet condition for this species on this surface"
            #If there is a Neumann condition on the potential, use that. Otherwise, the term depends on the solution for the potential
            if psurf in neumann.keys() and neumann[psurf][-1] is not None: #"if there IS a Neumann condition on the potential for this surface"
              expr=self.nbcs[(psurf,self.Nspecies)] #use the Neumann condition
            else:
              expr=fem.inner(self.n,fem.grad(Phi)) #depends on the solution
            #Store the expression for this term
            self.hbcs[(psurf,s)]=expr

    #Initial Conditions and Guess
    guesstup=self.species.initconc+[self.conditions.initial_potential]
    self.u.interpolate(fem.Constant(guesstup))
    self.u_k=fem.interpolate(fem.Constant(guesstup),self.V)
    u_klist=fem.split(self.u_k)
    c_klist=u_klist[:self.Nspecies]

    #Start time
    self.t=0.0
    #Calculate time step size
    self.dt=self.conditions.timedomain.stepsize

    #Weak Form
    beta=self.conditions.beta
    #Steady-state Nernst-Planck terms for each species
    #Volumetric terms
    weakforms=[]
    for s,c in enumerate(clist):
      if self.species.D[s] is not None:
        ##old weak form (based on Bin's code)
        ##J=-self.species.D[s]*(fem.grad(c)+self.conditions.beta*self.species.z[s]*c*fem.grad(Phi))
        ##wkform=fem.inner(J,fem.grad(vlist[s]))*fem.dx
        ##weakforms.append(wkform)
        term1=fem.inner(fem.grad(c),fem.grad(vlist[s]))*fem.dx
        if self.species.z[s] != 0:
          term2=beta*self.species.z[s]*c*fem.inner(fem.grad(vlist[s]),fem.grad(Phi))*fem.dx
          weakforms.append(-self.species.D[s]*(term1+term2))
        else:
          weakforms.append(-self.species.D[s]*term1)
    #Boundary terms for Neumann conditions
    for tup,expr in self.nbcs.items():
      psurf,i = tup
      if i < self.Nspecies:
        bterm = expr*vlist[i]*self.ds(psurf)
      else:
        bterm = -UN.eps_0*self.conditions.eps_r*expr*vPhi*self.ds(psurf)
      weakforms.append(bterm)
    #Hybrid boundary term
    for tup,expr in self.hbcs.items():
      psurf,s=tup
      if self.species.z[s] != 0:
        bterm = beta*self.species.z[s]*expr*clist[s]*vlist[s]*self.ds(psurf)
        weakforms.append(bterm)
    #Poisson terms
    poissonL=UN.eps_0*self.conditions.eps_r*fem.inner(fem.grad(Phi),fem.grad(vPhi))*fem.dx
    terms=[self.species.z[i]*clist[i] for i in range(self.species.N) if self.species.z[i] != 0]
    if len(terms)>0:
      poissonR=sum(terms)*vPhi*fem.dx
      #Add up to Steady-State PNP
      FF_ss=sum(weakforms)+poissonL-poissonR
    else:
      FF_ss=sum(weakforms)+poissonL
    #Time-dependent terms
    tdweaks=[]
    for i,c in enumerate(clist):
      term1=c*vlist[i]*fem.dx
      term2=c_klist[i]*vlist[i]*fem.dx
      tdweaks.append(term1-term2)
    #Reaction terms
    rxnweaks=[]
    for i,c in enumerate(clist):
      for j in range(self.reactions.N):
        if self.reactions.stoichio[j][i] != 0:
          termconst=self.reactions.stoichio[j][i]*self.reactions.constants[j]
          rxf=getattr(self,self.reactions.functions[j])
          term=termconst*rxf(*clist)*vlist[i]*fem.dx
          rxnweaks.append(term)
    #Put it all together
    if self.reactions.N > 0:
      self.FF=sum(tdweaks)-self.dt*FF_ss-self.dt*sum(rxnweaks)
    else:
      self.FF=sum(tdweaks)-self.dt*FF_ss
    #Take derivative
    self.J=fem.derivative(self.FF,self.u)

  def stopnow(self):
    """Check if stopping criterion is met.

    Returns True if time to stop, False otherwise"""

    criterion=self.conditions.timedomain.stopping
    if criterion.numsteps is not None:
      return self.k >= criterion.numsteps
    elif criterion.t_end is not None:
      return self.t >= criterion.t_end
    else:
      raise Exception("Unknown stopping criterion.")

  def run(self):
    """Do the time steps"""

    #Formulate problem and solver
    problem = fem.NonlinearVariationalProblem(self.FF, self.u, bcs=self.bcs, J=self.J)
    solver = fem.NonlinearVariationalSolver(problem)
    #Initialize time-domain output
    self.process_output_commands('datasteps')
    #Do the steps
    self.k=0
    while not self.stopnow():
      #Time at next step
      self.t+=self.dt
      for itm in self.timedeps:
        itm.t=self.t
      #Solve for next step
      solver.solve()
      #Output time-dependent values
      self.process_output_commands('datasteps')
      #Next step is now previous step
      self.u_k.assign(self.u)
      self.k+=1

simulatorclasses={'tdpnp_unhomog':TDPNPUSimulator}
