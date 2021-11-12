"""Create a function by interpolating from point data in a CSV file"""

#Standard library

#Site packages
import numpy as np
import fenics as fem
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.nested import WithNested
from . import simrequest
from . import meshinfo
from ..requesthandler import logging

logger=logging.getLogger(__name__)

_InterpolationConditions_props_schema_yaml="""#InterpolationConditions
functionname: {type: string}
coordcolumns:
  type: array
  items: {type: string}
valuecolumn: {type: string}
outcolumns:
  type: array
  items: {type: string}
boundaryvalue:
  anyOf:
    - {type: number}
    - {type: "null"}
nonboundary_fillvalue:
  anyOf:
    - {type: number}
    - {type: "null"}
"""

class InterpolationSimulator(simrequest.SimulationRequest):
  """Simulator for projecting an expression into a function space

  The point data must be loaded with a command in ``loaddata``,
  and saved in the attribute ``pointdata``.
  The result of the interpolation is saved in the ``soln`` attribute.
  
  Conditions:

    - functionname = optional string specifying name for projected function (not the attribute)
    - boundaryvalue = required function value to be used at mesh boundaries

      Use the value ``None`` to not define a boundary value.
      Note that this will result in NaN values unless the convex hull of the input points
      contains the entire mesh.

    - coordcolumns = optional list of names (as strings) for the columns containing X,Y,and Z, respectively.
    
      If there are only two dimensions, just leave off the column for Z.
      The default value is [x,y,z].

    - valuecolumn = optional name, as string, of the column containing the function value. (Defaults to 'f')

    - outcolumns = optional list of names (as strings) for the columns of the output dataframe.

    - nonboundary_fillvalue = optional "fill value" to specify for the interpolator even in the case when no boundary value is used

      If ``boundaryvalue`` is not ``None``, then this parameter is ignored.
      Defaults to ``nan``.

      Note that column names are all case-sensitive."""

  _validation_schema=simrequest.SimulationRequest.update_schema(
      _InterpolationConditions_props_schema_yaml,'properties.conditions.properties')
  _validation_schema.get_nested('properties.conditions.required').append('boundaryvalue')
    
  def run_sim(self):

    #For convenience
    self.conditions_processed=WithNested(**self.conditions)
    conditions=self.conditions_processed
    conditions.family=conditions.get_nested_default("family","P")
    d=self.meshinfo.spatial_dims()
    meshcoords=self.meshinfo.coordinates()
    Nmesh=meshcoords.shape[0]

    #Function space, and its coordinates
    self.V = fem.FunctionSpace(self.meshinfo.mesh,conditions.family, conditions.elementorder)
    Ndof_pts=self.V.dim()
    dofinfo=meshinfo.DOFInfo(self.meshinfo,self.V)
    dofcoords=dofinfo.coordinates()

    #Load the input data
    self.process_load_commands()
    df=self.pointdata

    #Get the function name
    functionname=getattr(conditions,'functionname','projected')

    #Get the column names for the input data
    coordcolumns=getattr(conditions,'coordcolumns',['x','y','z'])
    valuecolumn=getattr(conditions,'valuecolumn','f')
    outcolumns=getattr(conditions,'outcolumns',['dof_x','dof_y','dof_z','value'])

    #Input coordinates and function values, as separate arrays
    inpts=df.loc[:,coordcolumns].values
    invals=df.loc[:,[valuecolumn]].values.flatten()

    #Handle boundary value, if provided
    if conditions.boundaryvalue is None:

      data_pts=inpts
      data_vals=invals

      fillvalue=getattr(conditions,'nonboundary_fillvalue',np.nan)

    else:

      fillvalue = conditions.boundaryvalue

      #Get boundary and non-boundary DOF coordinates
      boundary_pts, nonbound_pts = dofinfo.boundary_points()

      #Array of boundary values
      boundary_vals=np.array([conditions.boundaryvalue]*boundary_pts.shape[0])

      #Combine boundary and input points
      data_pts=np.vstack([inpts,boundary_pts])
      data_vals=np.hstack([invals,boundary_vals])

    #Set up the interpolator
    logger.startTimer("interp_setup",request_name=getattr(self,"name",None))
    ilator=LinearNDInterpolator(data_pts,data_vals,fill_value=fillvalue,rescale=False)
    logger.stopTimer("interp_setup",request_name=getattr(self,"name",None))
    self.interp_setup_timer=logger.timers["interp_setup"]
    #Interpolate at all dof points
    logger.startTimer("interp_run",request_name=getattr(self,"name",None))
    dofvals=ilator(dofcoords)
    logger.stopTimer("interp_run",request_name=getattr(self,"name",None))
    self.interp_run_timer=logger.timers["interp_run"]
    #Store results for separate output if desired
    results_arr=np.hstack([dofcoords,np.reshape(dofvals,(dofvals.shape[0],1))])
    self.results=pd.DataFrame(results_arr,columns=outcolumns)
    #Define function from the interpolated dof values
    self.soln=fem.Function(self.V,name=functionname)
    junk=self.soln.vector().set_local(dofvals)
    del junk

    #Done
    return

  def compute_residual_errors(self,dfpath="pointdata",funcattr="soln",outattr="residuals"):
    df=self.get_nested(dfpath)
    func=self.get_nested(funcattr)
    #Get the column names for the input data
    coordcolumns=self.get_nested_default('conditions.coordcolumns',['x','y','z'])
    valuecolumn=self.get_nested_default('conditions.valuecolumn','f')
    #Set up the output dataframe
    resid=df.copy()
    #Compute the function value at each input point
    def calcvalue(row):
      pt=[row[c] for c in coordcolumns]
      try:
        res=func(*pt)
      except RuntimeError:
        res=np.nan #point is not inside mesh; use nan
      return res
    resid['interpolated']=df.apply(calcvalue,axis=1)
    #Compute the residuals
    def calcresidual(row):
      return row['interpolated']-row[valuecolumn]
    resid['error']=resid.apply(calcresidual,axis=1)
    #Store result
    self.set_nested(outattr,resid)
    return

#Register for loading from yaml
yaml_manager.register_classes([InterpolationSimulator])
