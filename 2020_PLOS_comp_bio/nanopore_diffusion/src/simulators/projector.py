"""Project an expression into a function space and save the result."""

#Standard library

#Site packages
import fenics as fem

#Local
import simulator_general

class Conditions(simulator_general.GenericConditions):
  """Condition defnitions for use with ProjectingSimulator

  Attributes:
  
    - functiontype = optional string specifying type of function:
    
      - 'scalar' (default) for a scalar function
      - 'vector' for a vector function
      
    - expression = required string, the expression to project
    
      Note that the python ``format`` method is called on the string, using the mesh metadata as the keyword arguments.
      This allows the expression to reference variables defining the mesh structure, without using FEniCS parameters.
      
    - parameters = optional list of parameters names to use in the expression, empty for no parameters.
      
      Note that the values of these parameters are taken from the Simulator's ``info`` attribute."""
  __slots__=['functiontype','expression']

class ProjectingSimulator(simulator_general.GenericSimulator):
  """Simulator for Unhomogenized Fickian Diffusion

  Additional attributes not inherited from GenericSimulator:
  
    - V = function space
    - expr = FEniCS Expression object for the expression to project
    - soln = Function resulting from the projection

"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(ProjectingSimulator, self).__init__(modelparams)

    #Get conditions
    self.conditions=Conditions(**modelparams.conditions)

    #Requested Function space
    ftype=getattr(self.conditions,'functiontype','scalar').lower()
    if ftype == 'scalar':
      self.V = fem.FunctionSpace(self.meshinfo.mesh,'P',self.conditions.elementorder)
    elif ftype == 'vector':
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', self.conditions.elementorder)
    else:
      raise Exception("Invalid functiontype: %s"%ftype)

    #Apply mesh metadata to the expressions
    exprstr=self.conditions.expression.format(**self.meshinfo.metadata)
    
    #Create the parameters dictionary
    params={}
    for k in getattr(self.conditions,'parameters',[]):
      params[k]=self.info[k]
    
    #Create the expression object
    self.expr=fem.Expression(exprstr,element=self.V.ufl_element(),**params)

  def run(self):
    "Perform the projection"
    self.soln=fem.project(self.expr,self.V)
    return

simulatorclasses={'projector':ProjectingSimulator}
