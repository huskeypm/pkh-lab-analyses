
#Output functions originally for 3D pore problems

#Standard library
import os.path as osp

#Site packages
import fenics as fem
import numpy as np

#Local
import folderstructure as FS
import plotdata

def domain_volume(self,pcell=None,attrname='domain_volume'):
  """Compute volume of a subdomain"""
  if pcell is None:
    this_dx=self.dx
  else:
    this_dx=self.dx(pcell)
  self.results[attrname]=fem.assemble(fem.Constant(1)*this_dx)
  return

#TODO: change this to store the flux in a specified attribute, use another function to save to VTK file
def fluxfield(self,filename, solnattr='soln', idx=None, fluxattr='flux', D_bulk=None):
  """Flux as vector field (new attribute, and VTK file)

  TODO: update arguments information

  Arguments:
    filename = name of output file, as string
      File will be created in the output directory (self.outdir)
  Required attributes:
    conditions.D_bulk = bulk diffusion constant for the medium
    outdir = output directory, as string
    soln = FEniCS Function containing solution
    V_vec = FEniCS VectorFunctionSpace
  New attributes:
    flux = flux field as a MeshFunction
  No return value.
  Output file is written."""
  if D_bulk is None:
    D_bulk=self.conditions.D_bulk
  soln=getattr(self,solnattr)
  if idx is not None:
    soln=soln[idx]
  expr=fem.Constant(-D_bulk)*fem.grad(soln)
  fluxres=fem.project(expr,self.V_vec,solver_type="cg",preconditioner_type="amg") #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)
  setattr(self,fluxattr,fluxres)
  if filename is not None:
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

def effective_diffusion(self,name,totflux_name):
  """Calculate effective diffusion constant
  Arguments:
    name = name for storage in the results dictionary
    totflux_name = name of previously calculated total flux in results dictionary
      This requires a previous call to fluxintegral.
  Required attributes:
    results[toflux_name] = result from previous call to fluxintegral
  No new attributes.
  New item added to results dictionary.
  No return value.
  No output files.
  CAUTION: this method is not generalized. It only works for certain geometries.
  In fact, probably only one: the body-centered nanopore.
  For the face-centered one, you'd need to sum the results of integration over two facets."""
  quarter_area = self.meshinfo.metadata['cell_area']/4
  zvals=[self.meshinfo.metadata['Z2'], self.meshinfo.metadata['Z3']]
  samples=[self.soln(0,0,zv) for zv in zvals]
  delta=samples[1]-samples[0]
  delta_z=zvals[1]-zvals[0]
  Deff=float(self.results[totflux_name]/quarter_area*delta_z/delta)
  self.results[name]=Deff
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

def volfrac(self,name):
  """Calculate free volume fraction
  Arguments:
    name = name for storage in the results dictionary
  Required attributes:
    mesh_metadata = dictionary of mesh metadata
      MUST contain 'pore_area' and 'cell_area' keys
  No new attributes.
  New item added to results dictionary.
  No return value.
  No output files."""
  self.results[name]=self.meshinfo.metadata['pore_area']/self.meshinfo.metadata['cell_area']
  return

#We could generalize this by specifying:
# - a center point
# - an axis of rotation defining the plane
# - a point to define the orientation of theta=0
def profile_radial(self,zval,theta,num,plotname,label,attrname='soln',idx=None):
  """Data for plot of solution along radial line at specified Z, in specified direction
  Arguments:
    zval = 
    theta = theta-angle in degrees from x-axis, as float
    num = number of sampled points
    indep = index of the coordinate parameter to use as the independent variable for the plot (zero-based)
    plotname = name of plot in outdata.plots, as string
    label = series label to assign, as string
    attrname = name of attribute to output, as string, defaults to 'soln'
  Required attributes:
    outdata = instance of OutData
    mesh_metadata = dictionary of mesh metadata
      MUST contain radius under key 'R'
  No new attributes.
  Nothing added to results dictionary.
  No return value.
  Series is added to `outdata.plots`."""
  #Calculate start and end locations
  startloc=(0,0,zval)
  rads=np.radians(theta)
  xend=self.meshinfo.metadata['R']*np.cos(rads)
  yend=self.meshinfo.metadata['R']*np.sin(rads)
  endloc=(xend,yend,zval)
  #Call line_profile to finish
  self.line_profile(startloc,endloc,num,plotname,label,attrname,None,idx)
