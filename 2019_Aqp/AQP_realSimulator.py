
from __future__ import division
from fenics import *
import numpy as np
import ufl
from scipy.optimize import fsolve

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import cPickle as pickle


def wrapf(f):
# takes either Function/uflExpression or float/tuple and wraps with Constant in the latter case
# for easy specification of Dirichlet or Neumann data
    if isinstance(f, GenericFunction) or isinstance(f, ufl.core.expr.Expr):
        return f
    elif isinstance(f, float) or isinstance(f, int) or isinstance(f, tuple):
        return Constant(f)
    else:
        dolfin_error(__name__+".py", "use given function",
            "Function is of unexpected type '%s'." % (type(f),))

def _invert_dict(d):
    if isinstance(d,dict):
        return {i:s for s,t in d.items() for i in t}




Dk = 1.96e-9
Dcl = 2.03e-9 # m^2/s  --Cl-
zh=1
zoh=-1
Dh=9.31e-9 # ---H+
Doh=5.3e-9 # ---OH-

  #---------------------------
  # Basic constants
  #---------------------------
F = 96485.3365 # C/mol
eps0 = 8.854187817e-12 # C/V-m
epsr = 78.5
kb = 1.380648e-23
T = 298 #K
charge = 1.60217e-19 #C --- electron charge

pH=7.4
nm=1e-9
rr = 2 # radius of pore
dd = 1 # thickness of buffering shell
#Btot=100 # ratio of [buffer]/[CF]
cKCl = 1.
meshfile ="/home/AD/srbl226/aqui/TEST.hdf5"
def runPNP(
  meshfile = meshfile,
  num = 4,
  bR = 3,
  cKCl = cKCl,
  phiPore = -10.
  ):
  bR=bR
  #aqD = aqD*nm
  phiAQ=0
  phiPore=0
  #print "phiAQ is: ", phiAQ
  hdf5Name = meshfile
  cCl = cKCl
  cK = cKCl
  #phi0 = phi
  #length = 10*nm
  #radius = rr*nm
  #spacing = space *nm
  #thickness = dd*nm

  #--------------------------------------------------

  zk = 0#1
  zcl = 0#-1 # Cl

  ch0 = 10**(-pH+3) #mol/m^3  initial [H]
  coh0 = 10**(-11+pH) #mol/m^3 inital [OH]
  ccl0 = cCl
  ck0 = cK
  parameters["allow_extrapolation"] = True
  parameters["ghost_mode"] = "shared_facet"
  mesh = Mesh()
  facets = MeshFunction("size_t",mesh)
  cells = MeshFunction("size_t",mesh)
  hdf5=HDF5File(mesh.mpi_comm(),hdf5Name,'r')
  hdf5.read(mesh,'mesh',False)
  hdf5.read(facets,'facets')
  hdf5.read(cells,'cells')
  hdf5.close()

  ds = Measure('exterior_facet',subdomain_data=facets, domain=mesh) # Interior surface integration
  dS = Measure('interior_facet',subdomain_data=facets, domain=mesh)
  dx = Measure("dx",subdomain_data=cells,domain=mesh)

  #Define MixedFunctionSpace--
  #--- for dolfin 1.6.0 (grimm)
  #P = FunctionSpace(mesh,"P",1)
  #V = MixedFunctionSpace([P,P,P,P,P,P])
  ele = FiniteElement('CG',mesh.ufl_cell(),1)
  mele = MixedElement([ele,ele,ele,ele,ele])
  V = FunctionSpace(mesh,mele)


  u = Function(V)
  u.interpolate(Constant((ck0,ccl0,ch0,coh0,0.0)))

  ck, ccl, ch, coh, v = split(u)
  ckm, cclp, chp, cohp, vv = TestFunctions(V)

  #-----------------------------------
  # Write the variational forms
  #----------------------------------
  # flux of ions
  Jm = -Dk*(grad(ck) + charge*zk/(kb*T)*ck*grad(v))
  Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
  Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
  Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))

  # get the LHS
  #  -div(J) = 0 ==> J*grad(testf)*dx = 0
  # here, we consider dx(0) is unbuffered region and dx(1) is bufferd region
  aJm = inner(Jm,grad(ckm))*dx
  aJp = inner(Jp,grad(cclp))*dx
  aJh = inner(Jh,grad(chp))*dx
  aJoh = inner(Joh,grad(cohp))*dx

  #--LHS and RHS of poisson equation
  aPoissonL = inner(grad(v),grad(vv))*dx
  aPoissonR = F/(eps0*epsr)*(ck*zk + ccl*zcl+ch*zh +coh*zoh)*vv*dx


  FF = aJm + aJp + aJh + aJoh  + aPoissonL - aPoissonR
  J = derivative(FF, u)



  #--------Boundary Conditions--------------------------
  #-- Ground Potential at the two ends of reservoirs
  bc1 = DirichletBC(V.sub(4),0,facets,3)
  bc2 = DirichletBC(V.sub(4),0,facets,4)

  #---------------------------------------
  # assigin boundary condition for K+ and Cl-
  bc3 = DirichletBC(V.sub(0),0,facets,4) #----> Assign a 0 [K+] at the botton reservor
  bc4 = DirichletBC(V.sub(0),ck0,facets,3)
  bc5 = DirichletBC(V.sub(4),Constant(0),facets,6)
  bc5x = DirichletBC(V.sub(1),0,facets,4)
  bc6 = DirichletBC(V.sub(4),Constant(0),facets,5)
  bc6x = DirichletBC(V.sub(1),ccl0,facets,3)
  # assign boundary condition for H+ and OH-
  bc7 = DirichletBC(V.sub(2),ch0,facets,4)
  bc8 = DirichletBC(V.sub(2),ch0,facets,3)
  bc9 = DirichletBC(V.sub(3),coh0,facets,4)
  bc10 = DirichletBC(V.sub(3),coh0,facets,3)
  bcx = DirichletBC(V.sub(4),Constant(0),facets,2)


  #-------------------------------------
  bcc = [bc1,bc2,bc3,bc4,bc5,bc5x,bc6,bc6x,bc7,bc8,bc9,bc10,bcx]
  #-------------------
  # Solve the problem
  #--------------------

  problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J) #deleted FF preceding u
  solver = NonlinearVariationalSolver(problem)



  (iter, converged) = solver.solve()
  #print "This one solved \n\n\n\n\n\n\n\n\n\n !!!!!!!!!!!!!"

  ck_u,ccl_u,ch_u,coh_u,v_u = u.split(True)



  TT = ck_u.function_space()
  degree = TT.ufl_element().degree()
  W = VectorFunctionSpace(mesh,'P',degree)
  #print "step before fluxck"
  fluxck = project(grad(ck_u)*Constant(Dk),W)
  #print "step after fluxck"

  v1file = File("solutions/AQ_{}_{}.pvd".format(cK,num))
  v1file << ck_u
  #v2file = File("solutions/ckflux_%s.pvd"%(str(cK)))
  #v2file << fluxck

  import os
  myPath = os.path.abspath(__file__)
  path = os.path.abspath(os.path.join(myPath,'..'))
  path = path+"/"
  hdf5=HDF5File(mesh.mpi_comm(),path+"AQ_{}_solution.hdf5".format(num),'w')
  hdf5.write(ck_u,"solution")
  hdf5.close()


  n = FacetNormal(mesh)
  area1 = assemble(Constant(1.0)*ds(1))#,domain=facets))
  #print "area1", area1
  area2 = assemble(Constant(1.0)*ds(2))#,domain=facets))
  #print "area2", area2
  area3 = assemble(Constant(1.0)*ds(3))#,domain=facets))
  #print "area3", area3
  area4 = assemble(Constant(1.0)*ds(4))#,domain=facets))
  #print "area4", area4
  area5 = assemble(Constant(1.0)*ds(5))#,domain=facets))
  #print "area5", area5
  area6 = assemble(Constant(1.0)*ds(6))#,domain=facets))
  #print "area6", area6
  area7 = assemble(Constant(1.0)*dS(7))#,domain=facets))
  area8 = assemble(Constant(1.0)*dS(8))#,domain=facets))
  #print "area7 ", area7
  #print "area8 ",area8
  #print "not broken yet"
  #print "pre dot"


  fl = dot(fluxck,n)
  flux_bot = assemble(fl('-')*dS(7))
  flux_top = assemble(fl('+')*dS(8))
  avgf = (flux_top)/(area1)
  avgbot = (flux_bot)/(area1)
  length =1e-8
  RevH = 1e-8
  Deff = avgf*(length + 2*RevH)/ck0#/Dk
  Deffbot = avgbot*(length + 2*RevH)/ck0
  #print "length", length
  #print "RevH", RevH
  #print "Deff", Deff
  file = open("AQ_{}.txt".format(num), "w")
  file.write("Deff {:3}\n".format(Deff))
  file.write("Deff {:3}\n".format(Deffbot))
  #file.write("{:3}".format(Deff))
  file.write("VolFrac {:3}".format(area5/area3))
  #file.write("{:3}".format(area7/area3))
  file.close()
  return Deff

###############
# Commenting out hte above until I get sigmaS working
  #pickle.dump(Results, open("Ca_%s_%s_%s.p"%(str(R),str(length/nm),str(cCaCl2)),"wb"))

#!/usr/bin/env python
import sys
#
# Message printed when program run without arguments
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose:
  Bin's solver

Usage:
"""
  msg+="mpirun -np #proc python  %s -validation" % (scriptName)
  msg+="""


Notes:
  Pay attention to the units in the input arguments!

"""
  return msg

#
# MAIN routine executed when launching this script from command line
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  1
  #  #print "arg"

  # Loops over each argument in the command line
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      runPNP()
      quit()
    if(arg=="-run"):
      arg1=sys.argv[i+1]
      arg2=np.float(sys.argv[i+2])
      arg3=np.float(sys.argv[i+3])
      arg4=np.float(sys.argv[i+4])
      runPNP(meshfile=arg1, phiPore=arg2,aqD=arg3,cKCl=arg4) #same as at the pickle line, justing doing radii
      quit()
    if(arg=="-runner"):
      arg1=sys.argv[i+1]
      arg2=int(sys.argv[i+2])
      #arg3=float(sys.argv[i+3])
      arg3=1#int(sys.argv[i+4])
      runPNP(meshfile=arg1,num=arg2, cKCl = arg3) #value provided is in mM = mol/m^3
      quit()
