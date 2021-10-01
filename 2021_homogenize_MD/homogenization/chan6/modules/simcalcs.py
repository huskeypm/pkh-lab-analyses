"""Post-processing calculations"""

# import numpy as np
import fenics as fem

def calc_delta(self, vmin, vmax, outattr):
  """Compute a delta from min and max

  Arguments:
  
    - vmin = minimum value
    - vmax = maximum value
    - outattr = attribute path for storing result"""
  dv=self.get_stored(vmax)-self.get_stored(vmin)
  self.set_nested(outattr,dv)
  return

def calc_ratio(self, numerator, denominator, outattr):
  a=self.get_stored(numerator)
  b=self.get_stored(denominator)
  res=a/b
  self.set_nested(outattr, res)

def calc_full_thickness(self, width_nm, slab_thickness, outattr):
  """Compute the unit cell thickness including slabs

  Arguments:
  
    - width_nm = channel width in nm
    - slab_thickness = slab thickness in nm
    - outattr = attribute path for storing area result"""
  chan_width=self.get_stored(width_nm)
  t_slab=self.get_stored(slab_thickness)
  full_t=chan_width + 2*t_slab
  self.set_nested(outattr, full_t)
  return

def effective_D(self,outattr,fluxattr,areaattr,start_conc_attr,end_conc_attr,delta_s_attr):
  """Calculate effective diffusion constant, using concentration averages rather than point values

  Arguments:

    - outattr = attribute path for storage of result
    - fluxattr = attribute path to previously calculated total

        This requires a previous call to fluxintegral (or similar).

    - areaattr = attribute path to previously calculated area in results dictionary

        This requires a previous call to facet_area (or similar).

    - start_conc_attr = attribute path to concentration value at starting boundary
    - end_conc_attr = attribute path to concentration value at ending boundary
    - delta_s_attr = attribute path to value of Delta s for effective D calculation

  No return value.
  No output files."""
  #Get the values from attributes
  totflux=self.get_nested(fluxattr)
  area=self.get_nested(areaattr)
  startconc=self.get_nested(start_conc_attr)
  endconc=self.get_nested(end_conc_attr)
  delta_s=self.get_nested(delta_s_attr)
  #Calculate the change in concentration between the two points
  delta_c=endconc-startconc
  #Calculate diffusion constant
  Deff=float(totflux/area*delta_s/delta_c)
  #Store result
  self.set_nested(outattr,Deff)
  return


#List of functions to be bound as methods
request_methods=[calc_ratio, calc_delta, calc_full_thickness, effective_D]