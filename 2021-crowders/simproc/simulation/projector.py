"""Project an expression into a function space and save the result."""

#Standard library

#Site packages
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.nested import WithNested
from . import simrequest

_ProjectorConditions_props_schema_yaml="""#ProjectorConditions
functionname: {type: string}
functiontype: {type: string}
projection_kwargs: {type: string}
"""

class ProjectionSimulator(simrequest.SimulationRequest):
  """Simulator for projecting an expression into a function space

  The expression must be loaded with a command in ``loaddata``.
  The result of the projection is saved in the ``soln`` attribute.
  
  Conditions:

    - functionname = string specifying name for projected function (not the attribute)

    - functiontype = optional string specifying type of function:
    
      - 'scalar' (default) for a scalar function
      - 'vector' for a vector function
      - 'matrix' for a rank-2 tensor function
      
    - projection_kwargs = optional dictionary of keyword arguments to the ``project`` function.
    
      This is used, for example, to set the linear solver and preconditioner."""

  _validation_schema=simrequest.SimulationRequest.update_schema(_ProjectorConditions_props_schema_yaml,'properties.conditions.properties')
  _validation_schema.set_nested('properties.conditions.required',[])
    
  def run_sim(self):

    #For convenience
    self.conditions_processed=WithNested(**self.conditions)
    conditions=self.conditions_processed
    conditions.family=conditions.get_nested_default("family","P")
    spatial_dims=self.meshinfo.mesh.geometry().dim()

    #Requested Function space
    ftype=getattr(conditions,'functiontype','scalar').lower()
    if ftype == 'scalar':
      self.V = fem.FunctionSpace(self.meshinfo.mesh, conditions.family, conditions.elementorder)
    elif ftype == 'vector':
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, conditions.family, conditions.elementorder)
    elif ftype == 'matrix':
      self.V = fem.TensorFunctionSpace(self.meshinfo.mesh, conditins.family, conditions.elementorder, (spatial_dims,spatial_dims))
    else:
      raise Exception("Invalid functiontype: %s"%ftype)

    #Load the expression
    self.process_load_commands()

    #Get the keyword arguments for projection
    projection_kwargs=getattr(conditions,'projection_kwargs',{})

    #Get the function name
    functionname=getattr(conditions,'functionname','projected')

    #Do the projection
    self.soln=fem.project(self.expr,self.V,**projection_kwargs)

    #Set the function name
    self.soln.rename(functionname,'')

#Register for loading from yaml
yaml_manager.register_classes([ProjectionSimulator])
