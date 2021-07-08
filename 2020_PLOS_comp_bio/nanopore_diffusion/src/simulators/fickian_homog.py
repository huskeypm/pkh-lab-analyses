"""Solve the homogenized Fickian diffusion problem
and extract the data needed for post-processing efforts"""

#Standard library

#Site packages
import fenics as fem

#Local
import unitsystem as UN
import common
import simulator_general

class Conditions(simulator_general.GenericConditions):
  """Condition defnitions for use with HomogFickianSimulator
  
  Additional attributes:
  
    - boundaries = list of physical lines representing boundaries, for construction of boundary terms"""
  __slots__=('boundaries')

class PeriodicBoundary(fem.SubDomain):
  """SubDomain subclass for Periodic boundary condition"""
  def __init__(self,xlims,ylims):
    """Arguments:
    
      - xlims = pair of x-values: (xmin,xmax)
      - ylims = pair of y-values: (ymin,ymax)"""
    super(PeriodicBoundary, self).__init__()
    self.left,  self.right = xlims
    self.bottom,self.top   = ylims
    self.xspan = self.right-self.left
    self.yspan = self.top-self.bottom
  # Left boundary is "target domain" G
  def inside(self, x, on_boundary):
    # return True if on left or bottom boundary AND NOT on one of the two corners 
    #  (self.left, self.top) and (self.right, self.bottom)
    return bool((fem.near(x[0], self.left) or fem.near(x[1], self.bottom)) and 
                (not ((fem.near(x[0], self.left) and fem.near(x[1], self.top)) or 
                      (fem.near(x[0], self.right) and fem.near(x[1], self.bottom)))) and on_boundary)
  def map(self, x, y):
    if fem.near(x[0], self.right) and fem.near(x[1], self.top):
      y[0] = x[0] - self.xspan
      y[1] = x[1] - self.yspan
    elif fem.near(x[0], self.right):
      y[0] = x[0] - self.xspan
      y[1] = x[1]
    else:   # fem.near(x[1], self.top)
      y[0] = x[0]
      y[1] = x[1] - self.yspan

class HomogFickianSimulator(simulator_general.GenericSimulator):
  """Simulator for Homogenized Fickian Diffusion
  
  Isotropy of the input diffusion constant is assumed.

  Additional attributes not inherited from GenericSimulator:

    - meshinfo = instance of simulator_general.MeshInfo
    - conditions = instance of HomogFickianConditions
    - V = FEniCS FunctionSpace on the mesh for chi
    - scalar_V = FEniCS FunctionSpace on the mesh for the diffusion constant
    - soln = FEniCS Function to store the result for chi
    - D = FEniCS Function to store the spatially varying diffusion constant
    - n = FEniCS FacetNormal
    - ds = FEniCS Measure for exterior facets
    - dx = FEniCS Measure for cells
    - solver = FEniCS LinearVariationalSolver"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(HomogFickianSimulator, self).__init__(modelparams)

    #Get conditions
    self.conditions=Conditions(**modelparams.conditions)

    #Periodic boundary condition
    xkeys=[k for k in self.meshinfo.metadata.keys() if k[0].upper()=='X']
    ykeys=[k for k in self.meshinfo.metadata.keys() if k[0].upper()=='Y']
    xvals=[self.meshinfo.metadata[k] for k in xkeys]
    yvals=[self.meshinfo.metadata[k] for k in ykeys]
    xlims=(min(xvals),max(xvals))
    ylims=(min(yvals),max(yvals))
    pbc = PeriodicBoundary(xlims,ylims)

    #Function Spaces and Functions
    #Function spaces
    self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', self.conditions.elementorder, constrained_domain=pbc)
    self.scalar_V = fem.FunctionSpace(self.meshinfo.mesh, 'P', self.conditions.elementorder)
    #Trial and Test Functions
    chi=fem.TrialFunction(self.V)
    v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V)

    #Point Dirichlet Boundary Condition
    # def point_bc_boundary(x, on_boundary):
    #   """Function for point Dirichlet boundary condition, as required by FEniCS"""
    #   return fem.near(x[0],xlims[1]/2.0) and fem.near(x[1],ylims[1]/2.0)
    # vec0=fem.Constant((0,0))
    # pointbc=fem.DirichletBC(self.V,vec0,point_bc_boundary,method='pointwise') 
    # self.bcs=[pointbc]
    self.bcs=[]

    #Load diffusion constant as a Function
    self.D=fem.Function(self.scalar_V)
    self.process_load_commands()
    
    #The index objects
    i=fem.i
    j=fem.j

    #Measure and normal for external boundaries
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Equation term dictionary
    eqnterms=simulator_general.EquationTermDict(simulator_general.EquationTerm)

    #Bilinear boundary terms
    for psurf in self.conditions.boundaries:
      termname="bilinear_boundary_%d"%psurf
      ufl=self.D*self.n[i]*chi[j].dx(i)*v[j]*self.ds(psurf)
      eqnterms.add(termname,ufl,bilinear=True)

    #Bilinear body terms
    termname="bilinear_body"
    ufl=-self.D*chi[j].dx(i)*v[j].dx(i)*self.dx
    eqnterms.add(termname,ufl,bilinear=True)

    #Linear boundary terms
    for psurf in self.conditions.boundaries:
      termname="linear_boundary_%d"%psurf
      ufl=self.D*self.n[i]*v[i]*self.ds(psurf)
      eqnterms.add(termname,ufl,bilinear=False)

    #Linear body terms
    termname="linear_body"
    ufl=-self.D*v[i].dx(i)*self.dx
    eqnterms.add(termname,ufl,bilinear=False)

    #FEniCS Problem and Solver
    a=eqnterms.sumterms(bilinear=True)
    L=eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(a,L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

  def run(self):
    "Run this simulation."
    self.solver.solve()

  def macroscale_diffusion(self,name="D_macro",attrname="soln",usevolume=None):
    """Perform the integral to obtain the homogenized diffusion constant
    
    Isotropy of the input D is assumed, but the output D may be anisotropic or even non-diagonal.

    Arguments:

      - name = optional, name for storage in the results dictionary, as string
      - attrname = optional, name for the attribute storing the solution (the result for chi)
      - usevolume = optional, name for the self.info key storing the volume to be used.
          
          If not provided, or ``None``, the volume will be taken as the total model volume.
          If provided, the name must exist in self.info.

    Required attributes (other than those from simulator_general):
    
      - the attribute given by attrname
      - dx = FEniCS Measure for cells
      - D = FEniCS Function to store the spatially varying diffusion constant
    
    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    if usevolume is None:
      volume=fem.assemble(fem.Constant(1)*self.dx)
    else:
      volume=self.info[usevolume]
    kdelta = lambda i,j: 1 if i==j else 0 #Kronecker delta
    soln=getattr(self,attrname)
    gradchi=fem.grad(soln)
    matr=[]
    for ii in range(2):
      row=[]
      for jj in range(2):
        term1=kdelta(ii,jj)*fem.assemble(self.D*self.dx)
        term2=fem.assemble(self.D*gradchi[jj,ii]*self.dx)
        val=(term1-term2)/volume
        row.append(val)
      matr.append(row)
    self.results[name]=matr

simulatorclasses={'fickian_homog':HomogFickianSimulator}