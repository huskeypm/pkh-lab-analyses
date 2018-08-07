"""
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

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
##
## Revisions
## 150416 PKH Upgraded to align with dolfin 1.4 as well as generalize python paths
## 

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from dolfin import *
import numpy as np
import sys
sys.path.append("/home/huskeypm/sources/jay/")
import poissonboltzmannSpecial as pbs
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
  poissonboltzmann.py -linear/-nonlinear/-finite/-wholedomain/-validation
  
Notes: 
  Guaranteed to be wrong for right now!!

Author:
  the computational scientist formally known as pete 



"""
## interpolation has been validated 
# Verified that Expression is being appropriately interpolated to mesh (3D sphere, Gamer output, 1D line and analytical)


minr = 12.
maxr = 22.

class params:
  def __init__(self):

    ## Parameters that are generally fixed
    # Standard units: A, mV, kcal, M
  
    # ec / 4 pi eps = 14.3840 [V Angstoms] 
    # --> eco4pieps = 14.384e3 [mV Angstroms]
    self.eco4pieps = 14.384e3 # conversion factor [mV Angstroms]
    #eco4pieps = 14.384 # [V Ang]
    self.kT = 0.593  # energy [kcal/mol]
    self.beta = 1/self.kT
    self.ec = 8.854187817e-12 # electric constant [C^2/(Jm)] !!! this is
    # vacuum perm, not ec
    self.M_TO_ANG = 1e-10
    self.J_TO_KCAL = 0.000239005736
    self.ec = self.ec / (self.M_TO_ANG * self.J_TO_KCAL) # ec [C^2/kcal A]
    self.epsilonExterior = 80. # dielectric constant in exterior []
    self.lambdaB = 7. # Bjerrum length [A] (Note: using value in Wiki http://en.wikipedia.org/wiki/Bjerrum_length)

    # other units
    F=96485.3365 # Faradays constant [C/mol]
    R = 8.3143   # Gas const [J/mol K] 
    T = 298.     # Temp [K] 
    self.RT_o_F=(R*T)/F*1e3  # [mV]   
    #print self.RT_o_F
    self.F_o_RT=1/self.RT_o_F# [1/mV]   

  
    ## Domain-specific parameters 
    self.dim = 2 # grid dimensions 
    #dim = 2  # 2d mesh 
    self.molRad=12.5 # radius of internal boundary in sphere
    #a=10.0 # radius of internal boundary in sphere
    self.molMarker = 2 # marker for molecular boundary 
    self.domMarker = 3  # marker for domain boundary 
    self.epsError = 0.001  # epsilson for 'error'
    self.epsError = 3.000  # epsilson for 'error' in selecting boundary
    self.domRad = 5. * self.molRad # radius of domain [A] (kind of, since square)
    
    ## System-specific parameters
    self.zLig  = -1        # unit charge of diffusing species (particle) 
    self.zProt = -3.       # unit charge of solute (big obstacle) 
    self.ionC = .150 # ion conc [M]
    self.ionRad=2. # ion radius [A]

    ## solver modes
    self.mode = "linear"# linear, nonlinear, finitesize 
    #mode = "nonlinear"
    #mode = "finite"    


  # ion
  def update(self):
    self.center = np.zeros(self.dim)      
    self.kappa = 0.328 * np.sqrt(self.ionC) # inverse debye length for monovalent specieds [1/A]
    self.ikappa  = 1/self.kappa # Debye length [A], "Intermoplecular and surface forces, Israelachvili" 
    #self.Fz_o_RT=self.zLig/self.RT_o_F     # [1/mV] 
    


parms = params()

## expression
# WARNING: get answer from Holst to make sure 
# Assumes charge is centered at origin
# WARNING: not validated for 2d 
def DebyeHuckelExpr(dim=3):
  # e^{ka}/(1+ka) Z ec / (4 pi eps eps0) * e^{-kr}/r
  parms.update()
  prefac=parms.zProt * np.exp(parms.kappa*parms.molRad) *parms.eco4pieps 
  prefac/=parms.epsilonExterior*(1+parms.kappa*parms.molRad)
  # f is in units [mV]
  # validated using a=10.0, ionC=0.1, Z=5 
  if(dim==3):
    f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])",prefac=prefac,k=parms.kappa)
  elif(dim==2):
    f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]+x[1]*x[1]))/sqrt(x[0]*x[0]+x[1]*x[1])",prefac=prefac,k=parms.kappa)      
  elif(dim==1):
    f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]))/sqrt(x[0]*x[0])",prefac=prefac,k=parms.kappa)      
  
  return f

## expression
# WARNING: get answer from Holst to make sure 
# Assumes charge is centered at origin
# WARNING: not validated for 2d 
# From http://en.wikipedia.org/wiki/DLVO_theory
def DLVOExpr(dim=3):
  # e^{ka}/(1+ka) Z ec / (4 pi eps eps0) * e^{-kr}/r
  parms.update()
  prefac=parms.zLig**2 * parms.lambdaB * (np.exp(parms.kappa*parms.molRad) *parms.eco4pieps)**2
  print prefac
  prefac/=parms.epsilonExterior*((1+parms.kappa*parms.molRad)**2)
  print prefac
  quit()
  prefac = 1.

  print "WARNING: this is not validated (z^2 due to prot*lig?)!!"

  # f is in units [mV]
  # validated using a=10.0, ionC=0.1, Z=5 
  if(dim==3):
    f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])",prefac=prefac,k=parms.kappa)
  elif(dim==2):
    f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]+x[1]*x[1]))/sqrt(x[0]*x[0]+x[1]*x[1])",prefac=prefac,k=parms.kappa)      
  elif(dim==1):
    f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]))/sqrt(x[0]*x[0])",prefac=prefac,k=parms.kappa)      
  
  #f = Constant(0.)
  return f


# share parms 
pbs.parms = parms


##
## BCs
## 

# covers anything inside the box boundaries 
class crowderBoundary(SubDomain):
  def inside(self,x,on_boundary):
    #result = np.linalg.norm(x-parms.center) > (-parms.epsError+parms.domRad)
    resultx = x[0] > (parms.mins[0]+parms.epsError) and x[0] < (parms.maxs[0]-parms.epsError)
    resulty = x[1] > (parms.mins[1]+parms.epsError) and x[1] < (parms.maxs[1]-parms.epsError)
    
    within = resultx and resulty and on_boundary
    if within:
      result = True
    else: 
      result = False
    #if on_boundary and result==False:
    #if result:
    #  print x, result,on_boundary 
    return result

class boxBoundary(SubDomain):
  def inside(self,x,on_boundary):
    #result = np.linalg.norm(x-parms.center) > (-parms.epsError+parms.domRad)
    resultx = x[0] > (parms.mins[0]+parms.epsError) and x[0] < (parms.maxs[0]-parms.epsError)
    resulty = x[1] > (parms.mins[1]+parms.epsError) and x[1] < (parms.maxs[1]-parms.epsError)

    within = resultx and resulty 
    if within:
      result = False
    else:
      result = True  
    result = result and on_boundary 

    #if result:
    #  print x

    return result



class molecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-parms.center) < (parms.epsError+parms.molRad)
    result = result and on_boundary
#    if result:
#      print np.linalg.norm(x)
#      print x
    return result      

class domainBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-parms.center) > (-parms.epsError+parms.domRad)
    result = result and on_boundary
    #print x
    #print result
    return result      

##
## Solver 
## 

# Should return potential in units [mV] 
def SolvePoissonBoltzmann(mesh,meshType="dolfin",boundaryPotential="DH",outerZero=True,boundaryCondition="molecular",outName="pbsolution.pvd"):
  #
  parms.dim = np.shape(mesh.coordinates())[1]
  parms.center = np.zeros(parms.dim)      
  parms.mins = np.min(mesh.coordinates(),axis=0)
  parms.maxs = np.max(mesh.coordinates(),axis=0)
  parms.update()


  V = FunctionSpace(mesh, "Lagrange", 1)

  ## Define boundary conditions
  subdomains = MeshFunction("size_t",mesh,parms.dim-1)
  subdomains.set_all(0)
  bcs = []

  # define BC on molecular boundary 
  if boundaryCondition=="molecular":
    print "Assuming your sphere is centered at 000"
    innerboundary = molecularBoundary()
    outerboundary = domainBoundary()
  elif boundaryCondition=="crowders":
    print "Assuming everything within box dims are crowders"
    innerboundary = crowderBoundary()
    outerboundary = boxBoundary()

  innerboundary.mark(subdomains,parms.molMarker)
  # PKH: should maybe use sympy later 
  # Analytical solution for linearized PBE for atom of radius R and charge q
  # see Eq 5.1 [1]	
  #M. Holst, N. Baker, and F. Wang, JCC, v21, n15,pp1319-1342, 2000
  #q = parms.ec * parms.z
  #k = parms.kappa/np.sqrt(parms.epsilonExterior)
  #u0 = Expression("q/epsilon*R*(1-k*R/(1+k*a))",\
  #  q=q,k=k,epsilon=parms.epsilonExterior,\
  #  R=parms.molRad,a=parms.molRad)    
  if(boundaryPotential=="DH"): 
    f = DebyeHuckelExpr(parms.dim)
  else:
    f = Constant(boundaryPotential)
  bcs.append(DirichletBC(V, f, subdomains,parms.molMarker))

  # test 
  m = Function(V)
  #bc = DirichletBC(V,Constant(1.),boundary)
  #bc.apply(m.vector())
  m.vector()[:]=0
  bcs[0].apply(m.vector())
  File("boundaryPotential.pvd") << m


  # define BC on domainboundary (potential should be zero at boundary) 
  if(outerZero):
    print "WARNING: this is only correct for when potential actually decays to 0"
    outerboundary.mark(subdomains,parms.domMarker)
    f = Constant(0.)
    bcs.append(DirichletBC(V, f, subdomains,parms.domMarker))
  
    m = Function(V)
    bc = DirichletBC(V,Constant(1.),subdomains,parms.domMarker)
    bc.apply(m.vector())
    File("outerBoundary.pvd") << m

  
  return PBEngine(mesh,V,subdomains,bcs,meshType,outName=outName)

# Not very general PB solver 
def PBEngine(mesh,V,subdomains,bcs,meshType="dolfin",outName="pbsolution.pvd"):
  
  ## Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  ## define form functions 
  import ufl
  namespace = ufl.__dict__.copy()
  
  # exterior
  # LHS  
  if(meshType=="dolfin"): 
    form = -1*parms.epsilonExterior*inner(grad(u), grad(v))*dx(domain=mesh)
  else:
    form = -1*parms.epsilonExterior*inner(grad(u), grad(v))*dx(1)

  # eps*grad(u)grad(v) = kappa^2 uv
  if(parms.mode=="linear"):
    if(meshType=="dolfin"): 
      form += -1*parms.epsilonExterior*parms.kappa*parms.kappa *u*v*dx
  # eps*grad(u)grad(v) = kappa^2 sinh(u)*v
  elif(parms.mode=="nonlinear"):
    form = pbs.Nonlinear()

  # eq (7) notes.pdf
  elif(parms.mode=="finite"):
    form = pbs.Finite()

  else:
    raise RuntimeError("What did you just tell me to do???")
  
  ### END
  #raise RuntimeError("")
  #form += Constant(0.) * v * dx
  form += Constant(0.)*v*dx # Never understood why I have to do this. Otherwise get error about incomplete form  
  #print "Solving %s form of PBE" % parms.mode
  x = Function(V)
  solve(lhs(form)==rhs(form), x, bcs)


  File(outName) << x
  
  return (V,x)




# Mostly to check that interpolations are correct
def ValidateDebyeHuckel(meshFile):
  f = DebyeHuckelExpr()

  ## 3D mesh interpolation
  #mesh = Mesh("sphere_mesh.xml.gz")
# mesh = UnitSphere(50)
# mesh.coordinates()[:] = 50*mesh.coordinates()[:]
# V = FunctionSpace(mesh,"CG",1)
# phi = interpolate(f,V)
# #(gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,-255:255:500j]
# (gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,minr:maxr:100j]
# interp = griddata(mesh.coordinates(),phi.vector(),(gx,gy,gz))
# interp[np.isnan(interp)]=0
# interp = np.reshape(interp,100)
# gz = np.reshape(gz,100)
  
  ## GAMER 
  meshg = Mesh(meshFile)                
  Vg = FunctionSpace(meshg,"CG",1)
  phig = interpolate(f,Vg)
  #(gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,-255:255:500j]
  (gxg,gyg,gzg) = np.mgrid[0:0:1j,0:0:1j,minr:maxr:100j]
  interpg = griddata(meshg.coordinates(),phig.vector(),(gxg,gyg,gzg))
  interpg[np.isnan(interpg)]=0
  interpg = np.reshape(interpg,100)
  gzg = np.reshape(gzg,100)
  
  
  
  ## 1D mesh
  f = DebyeHuckelExpr(dim=1)
  mesh1 = UnitIntervalMesh(150) # 0..1
  m=mesh1.coordinates()
  mesh1.coordinates()[:]=m*(maxr-minr) + minr
  V1 = FunctionSpace(mesh1,"CG",1)
  phi1 = interpolate(f,V1)
  gz1= np.mgrid[minr:maxr:100j]
  interp1 = griddata(mesh1.coordinates(),phi1.vector(),gz1)
  
  ## analytical
  print "WARNING: not sure what Expression assumes if just one argument is passed in"
  def vec(r):
    return f(r)
  vecf = np.vectorize(vec)
  
  # double check
  # agrees with Expression 
  def func(r,prefac,k):
    x = np.array([r,0,0])
    p =prefac*np.exp(-k*np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    return p
  
  
  ## compare 
  plt.figure()
  #plt.plot(gz,interp,"r", label="3d UnitSphere")
  plt.plot(gzg,interpg,"r.", label="3d gamer")
  plt.plot(gz1,interp1,"g",label="1d UnitInterval")
  plt.plot(gz1,vecf(gz1),"b",label="analytical")
  plt.xlim([minr,maxr])
  #plt.legend(['a','b','c','d'])
  plt.legend(loc=2)                         
  title = "Debye-Huckel potential for %4.1f [A] sphere, 1/k= %4.1f [A], z = %d" %\
    (parms.molRad, parms.ikappa, parms.zProt)
  plt.title(title)
  plt.xlabel("r [A]")
  plt.ylabel("potential [mV]")
  plt.gcf().savefig("debyehuckel.png")

  
def ValidatePoissonBoltzmann(meshFile):

  ## GAMER 
  mesh = Mesh(meshFile)               
  (V,x) = SolvePoissonBoltzmann(mesh)

  ## validation of solution 
  # interpolate solution to line 
  (gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,minr:maxr:100j]
  interp = griddata(mesh.coordinates(),x.vector(),(gx,gy,gz))
  interp[np.isnan(interp)]=0
  interp = np.reshape(interp,100)
  gz = np.reshape(gz,100)
  
  # analytical 
  f = DebyeHuckelExpr()
  phi = interpolate(f,V)
  interpv = griddata(mesh.coordinates(),phi.vector(),(gx,gy,gz))
  interpv[np.isnan(interpv)]=0
  interpv = np.reshape(interpv,100)

  ## compare 
  plt.figure()
  plt.plot(gz,interp,"r", label="Numerical")
  plt.plot(gz,interpv,"b",label="analytical")
  #plt.plot(gz1,vecf(gz1),"b",label="analytical")
  plt.xlim([minr,maxr])
  #plt.legend(['a','b','c','d'])
  plt.legend(loc=2)                         
  title = "Debye-Huckel potential for %4.1f [A] sphere, 1/k= %4.1f [A], z = %d" %\
    (parms.molRad, parms.ikappa, parms.zProt)
  plt.title(title)
  plt.xlabel("r [A]")
  plt.ylabel("potential [mV]")
  plt.gcf().savefig("debyehuckel_solve.png")

# Performs validation against DebyeHuckel equation 
def Validations():
  parms.dim = 3
  parms.center = np.zeros(parms.dim)
  # 15 ensures all points at boundary r=12.5 are found   
  # 200 ensures all points at boundary r=225 are found   
#  parms.innerBound = 15
#  parms.outerBound = 220
  meshFile = "./sphere_mesh.xml.gz"


  # if these fail, check that correct Debye Huckel expression 
  # is used for grid 'dimension' 
  ValidateDebyeHuckel(meshFile)
  ValidatePoissonBoltzmann(meshFile)

# from Graham equation for 1:1 electrolyte
# Israelachvili pg 310
# sigma = 0.117 sinh(psi0/51.4) sqrt([NaCl])
# NaCl [M], psi0 [mV] 
# sigma: surface Charge density [C/m^2] 
# Compare w Tbl 14.1 of Israelachvili  
def Grahame(sigma,cNaCl):
  psi=np.arcsinh(sigma/(0.117*np.sqrt(cNaCl)))*51.4
  return psi


      



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
    print arg
    #if(arg=="-mytest"):
    #  Mytest()
    #  quit()
    if(arg=="-nonlinear"):
      parms.mode="nonlinear"
    if(arg=="-finite"):
      parms.mode="finite"
    if(arg=="-wholedomain"):
      mode = "wholedomain"
      fileIn= "./example/sphere/sphere_2d_entire.xml"
    if(arg=="-linear"):
      parms.mode = "linear"
      parms.molRad = 1.5
      parms.domRad = 5. 
      #fileIn= "./example/sphere/sphere2d.xml"
      fileIn= "./example/2d/volFrac_0.27.xml"
    if(arg=="-validation"):
      Validations()
      quit()
  quit()


  if(mode=="wholedomain"):
    pbs.domainBoundary = domainBoundary
    pbs.molecularBoundary = molecularBoundary
    pbs.doWholeDomainPB(fileIn)
  elif(mode=="default"): 
    parms.epsError = 1.000  # epsilson for 'error' in selecting boundary
    parms.molRad=12.5 # actual radius is 
    parms.domRad=15.
    mesh = Mesh(fileIn)
    SolvePoissonBoltzmann(mesh)


