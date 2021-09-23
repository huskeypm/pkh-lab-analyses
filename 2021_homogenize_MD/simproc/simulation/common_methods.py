"""Module holding methods that are used by more than one simulator subclass, but not enough to be in simrequest itself."""

import numpy as np
import fenics as fem

def calcflux(self, solnattr='soln', idx=None, attrpath='flux', Dattr='Dbar_proj', funcname="flux"):
  """Flux as vector field (new attribute)"""
  soln=self.get_nested(solnattr)
  Dvalue=self.get_nested(Dattr)
  if idx is not None:
    soln=soln[idx]
    Dvalue=Dvalue[idx]
  expr=-Dvalue*fem.grad(soln)
  fluxres=fem.project(expr,self.V_vec,solver_type="cg",preconditioner_type="amg") #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)
  fluxres.rename(funcname,"calculated flux")
  self.set_nested(attrpath,fluxres)
  return

def fluxintegral(self,pfacet,attrpath,internal=False,fluxsign=None,normalvar=None,fluxattr='flux',idx=None):
  """Flux integral over specified facet

  Arguments:
    pfacet = physical facet number for flux measurement
    attrpath = attribute path for storage of the integral
    internal = boolean, default False, True to use internal boundary, False for external
    fluxsign = '+' or '-' to specify which direction normal to the facet for flux calculation
      Required only if internal==True
    fluxattr = name of attribute containing the flux field for integration
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
  this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
  totflux=fem.assemble(fem.dot(getattr(self,fluxattr),this_n)*this_ds(pfacet))
  self.set_nested(attrpath,totflux)
  return

def effective_D(self,outattr,fluxattr,areaattr,startloc,endloc,solnattr='soln',idx=None):
  """Calculate effective diffusion constant

  Arguments:

    - outattr = attribute path for storage of result
    - fluxattr = attribute path to previously calculated total

        This requires a previous call to fluxintegral.

    - areaattr = attribute path to previously calculated area in results dictionary

        This requires a previous call to facet_area.

    - startloc = argument to get_pointcoords for start of line
    - endloc = argument to get_pointcoords for end of line
    - solnattr = optional, attribute path to the concentration solution
    - idx = index of the solution field to write out, None (default) if not a sequence

  No return value.
  No output files."""
  #Get the flux and area values
  totflux=self.get_nested(fluxattr)
  area=self.get_nested(areaattr)
  #Get the object with the solution data
  vals=self.get_nested(solnattr)
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
  Deff=float(totflux/area*delta_s/delta_c)
  #Store result
  self.set_nested(outattr,Deff)
  return
