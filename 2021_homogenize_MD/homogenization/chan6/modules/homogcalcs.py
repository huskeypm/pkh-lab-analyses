"""Post-processing calculations for homogenization simulations"""

import fenics as fem

def calc_ratio(self, numerator, denominator, outattr, use_stored=True):
  if use_stored:
    a=self.get_stored(numerator)
    b=self.get_stored(denominator)
  else:
    a=self.get_nested(numerator)
    b=self.get_nested(denominator)
  res=a/b
  self.set_nested(outattr, res)

def calc_product(self,factor1,factor2,outattr):
  a=self.get_stored(factor1)
  b=self.get_stored(factor2)
  res=a*b
  self.set_nested(outattr,res)

def calc_delta(self, vmin, vmax, outattr):
  """Compute a delta from min and max

  Arguments:
  
    - vmin = minimum value
    - vmax = maximum value
    - outattr = attribute path for storing result"""
  dv=self.get_stored(vmax)-self.get_stored(vmin)
  self.set_nested(outattr,dv)
  return

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

def project_exp_pot(self,outattr="exp_pot",funcname="exp_pot",solver="cg",precond="amg",):
  """Project the exponential of the potential, for the calculation of Xi

  Arguments:

    - outattr: attribute path for storing the result
    - funcname: name to be given to the projected function
    - solver: linear solver to be used for the projection operation
    - precond: preconditioner to be used for the projection operation

  """
  beta = self.conditions['beta']
  expr = fem.exp(-beta*self.potential)
  res = fem.project(expr,self.scalar_V,solver_type=solver,preconditioner_type=precond)
  res.rename(funcname,funcname)
  self.set_nested(outattr,res)

#List of functions to be bound as methods
request_methods=[calc_ratio, calc_product, calc_delta, calc_full_thickness, project_exp_pot]
