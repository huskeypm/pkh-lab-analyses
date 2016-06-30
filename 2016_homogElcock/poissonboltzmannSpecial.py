"""
----------------------------------
# ADD OWNERSHIP

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"""
from dolfin import *
import numpy as np
#import modelparameters

#
# Revisions
# Todo
# -validate all aspects of code 
#  The sinh expression initially comes out as zero, so it never looks like it is used in the form expression 
# -implement log function 
# -wholeDom prob. I think is behaving reasonably (need to visualize the potential in the exterior domain at a different dynamic range [0-10] than the interior domain)


msg="""
Purpose:
  Tempate for solving (non)-linear poisson boltzmann equation 

Usage:
  poissonboltzmann.py run 
  
Notes: 
  Guaranteed to be wrong for right now!!

Author:
  the computational scientist formally known as pete 



"""


print "WARNING: This code is no-where close to being validated (use at your own risk)"

# exterior domain 
class OmegaOutside(SubDomain):
    def inside(self, x, on_boundary):
      result = np.linalg.norm(x-params.center) > (-params.epsError+params.molRad)
      #print "outside"
      #print result 
      return result 

# interior domain 
class OmegaInside(SubDomain):
    def inside(self, x, on_boundary):
      result = np.linalg.norm(x-params.center) <=(-params.epsError+params.molRad)
      #print x
      #print params.center
      #print "inside"
      #print result 
      return result 

## special functions 
def sinh(x=0):
 print "Craptastic implementation. Need to use new version of fenics w sinh/Johans sympy model module"  
 return (exp(x)-exp(-x))/2.

def cosh(x=0):
 print "Craptastic implementation. Need to use new version of fenics w sinh/Johans sympy model module"  
 return (exp(x)+exp(-x))/2.


def doWholeDomainPB(fileIn):
  params.domRad = 10.  # radius of domain [A] (kind of, since square)
  scale = 2*params.domRad 
  ## mesh/functions 
  square=1 
  if square:
    mesh = UnitSquare(50,50)   
    mesh.coordinates()[:] = scale * mesh.coordinates()
    print "I lied, using hardcoded geom"
  else:
    # domain assigment doesn't seem to work correctly (fails on 'choose' funtion 
    mesh = Mesh("sphere_2d_entire.xml") 

  dim = 2
  V = FunctionSpace(mesh, "DG", 0)
  eps = Function(V)
  kappa= Function(V)


  ## point source 
  print "Point source is not functional"
  params.center = scale*np.array([0.5,0.5])
  p = Point(params.center[0],params.center[1])  
  #c = PointSource(V,p,magnitude=1) 
  q1 = 1. # [e]
  c = PointSource(V,p,q1)                 

  ## boundary/domain conditions 
  subdomains = MeshFunction('uint', mesh, dim)
  subdomainBoundary = MeshFunction('uint', mesh, dim-1)
  bcs = []
  boundary = domainBoundary()
  boundary.mark(subdomainBoundary,params.domMarker)
  u0 = Constant(0.0)
  bcs.append(DirichletBC(V, u0, subdomainBoundary,params.domMarker))

  ## discontinuous epsilon 
  # http://fenicsproject.org/documentation/tutorial/materials.html#working-with-two-subdomains
  subdomain0 = OmegaOutside()
  subdomain0.mark(subdomains, 0)
  subdomain1 = OmegaInside()
  subdomain1.mark(subdomains, 1)
 
  eps_values = [params.epsilonExterior, 1.]  # values of k in the two subdomains
  kappa_values = [params.epsilonExterior*params.kappa**2, params.epsError]  # values of k in the two subdomains
  #for cell_no in range(len(subdomains.array())):
  #  subdomain_no = subdomains.array()[cell_no]
  # eps.vector()[cell_no] = eps_values[subdomain_no]
  help = np.asarray(subdomains.array(), dtype=np.int32)
  eps.vector()[:] = np.choose(help, eps_values)

  kappa.vector()[:] = np.choose(help, kappa_values)
  #xkappa = eps


  ## Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  form = -1*eps*inner(grad(u), grad(v))*dx
#  form += eps *v*dx

  # eps*grad(u)grad(v) = kappa^2 uv
  form += -1*kappa *u*v*dx

  d= 1.0
  factor = 4 * np.pi * params.ec*params.ec / (params.kT)
  norm = 1/(d * np.sqrt(np.pi*2))
  form += Expression("factor*norm *exp(-( (x[0]-xC)*(x[0]-xC) + (x[1]-yC)*(x[1]-yC))/(2*d*d))",\
                      xC=params.center[0],yC=params.center[1],d=d,norm=norm,factor=factor)*v*dx
  A = lhs(form)
  b = rhs(form)

  # apply point charge 
   # not correct 
  #c.apply(b)
  

  x = Function(V)
  solve(A==b, x, bcs)
  # F = a
  # problem = NonlinearVariationalProblem(F, x, bcs=bc, J=J)
  # solver = NonlinearVariationalSolver(problem)
  #solver.parameters["linear_solver"] = "gmres"
  #solver.parameters["preconditioner"] = "ilu"
  #solver.solve()


  File("out.pvd") << x
  

def Nonlinear():
    print "Idiot, your weak form is non-linear in u - need to write newton iter!"
    print "Or check out pb work in smolhomog"
    quit()
    print "Using nonlinear code"
    namespace["sinh"] = sinh
    sinh_expr = eval(str(-1*params.epsilonExterior*params.kappa*params.kappa * sinh()), namespace, {"x":u})
    print "sinh_expr %s" % sinh_expr

    print "WARNING: not sure if ufl exp is reflected in the form (need to compare linear w non-linear"
    form += sinh_expr*v*dx
    print form 

    # Might need to cmopile this as UFL code, see pg 196 of dolfin manual
    # https://answers.launchpad.net/dolfin/+question/199782
    # WRONG form += -1*params.kappa*params.kappa *Expression("sinh(u)",u=u)*v*dx
    # WRONG form += -1*params.kappa*params.kappa *f(u)*v*dx

  # eq (7) notes.pdf
def Finite():
#  elif(params.mode=="finite"):
    print "Idiot, your weak form is non-linear in u - need to write newton iter!"
    quit()
    print "WARNING: hyperbolic sine not recognized"
    arg = params.z*params.beta*params.ec*u

    sinh_expr = eval(str(sinh()), namespace, {"x":arg})
    print "sinh_expr %s" % sinh_expr
    cosh_expr = eval(str(cosh()), namespace, {"x":arg})
    print "cosh_expr %s" % cosh_expr

    a= params.ionRad
    phi0 = 2*a*a*a*params.ionC     
    rhs1 = 8 *np.pi * params.z*params.ec*params.ionC
    rhs1 *= 1/params.epsilonExterior
    num = sinh_expr      
    denom = (1-phi0+phi0*cosh_expr)      
    rhs1 *= num/denom
    form += -1*rhs1*v*dx
    print form 

#sphere
if __name__ == "__main__":

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= "./example/sphere/sphere2d.xml"
  mode = "outer"
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-nonlinear"):
      params.mode="nonlinear"
      print "not functional"
      #quit()
    if(arg=="-wholedom"):
      mode = "wholedom"
      fileIn= "./example/sphere/sphere_2d_entire.xml"


   



  if(mode=="wholedom"):
    doWholeDomainPB(fileIn)
  else: 
    doOuterDomainPB(fileIn)


