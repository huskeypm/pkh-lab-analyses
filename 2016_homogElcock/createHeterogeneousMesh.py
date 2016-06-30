

## Load everying in 
import sys
root="/home/pmke226/sources/"
root="/home/AD/pmke226/sources/"
sys.path.append(root+"/smolhomog/")
sys.path.append(root+"/modified-pb/example/")
sys.path.append(root+"/cytoplasm/")              
import poissonboltzmann as pb 
import runLattice as rL
import homoglight as hl

from dolfin import * 
import numpy as np

##
## Assign parameeters 
##
parms = pb.parms
parms.kappa = 0.01 # 1/[A] 
length=1000. # Box length [A] (note: this might be overwritten in DefineCases)
eps = length/100.


class empty:pass

# Representative case
def DefineCases_Exp1():
  cases = dict()

  length = 1000.  # [A] 
  alllocs = np.array([[250,250],[750,750],[750,250],[250,750]])
  #cases["locs"] = locs
  allrads = np.ones( np.shape(alllocs)[0]) * 100.
  #cases["rads"] = rads
  print alllocs
  quit()
  
  # print postively charged points 
  points = empty()
  cases['positive'] = points
  #points.locs = [locs[0]]
  #points.rads = [rads[0]]
  inds = [0,1]
  points.alllocs = alllocs
  points.allrads = allrads
  points.inds = inds
  points.marker = 1
  points.V0 = 25.6  # mV
  
  points = empty()
  cases['negative'] = points
  #points.locs = [locs[1]]
  #points.rads = [rads[1]]
  inds = [2,3]
  points.alllocs = alllocs
  points.allrads = allrads
  points.inds = inds
  points.marker = 2
  points.V0 = -25.6  # mV
  

  return cases 

# Calculates effective diffusion constants for
# a system with crowders of mixed potentials 
def CalcDiffusivityMixedChargeMesh(
    cases,
    fileName = "temp_mesh",
    makeMesh=True,
    res=3., # [Mesh resolution] 
  ): 
  print "WARNING: seems to be ignoring my res= values"
  xmlName=fileName+".xml"  
  
  # runs homog on a lattice with non-uniformly distributed crowder locations
  # and radii 
  boxMin = [0,0]
  boxMax = [length,length] 
  
  if makeMesh:
    akey = cases.keys()[0]
    dxr,dyr,results = rL.HomogLattice(cases[akey].alllocs,
        res = res,      
        rads=cases[akey].allrads,          
        computeMargin=False,
        boxMin=boxMin,boxMax=boxMax,
        gmshName = fileName+".geo"
      )
  
  # Boundary 1 (near centroid) 
  class Crowder(SubDomain):
    def inside(self,x,on_boundary):
      result = np.linalg.norm( self.centroid - x ) < (self.radius + eps)
      #result  = x[0] < 0.5 + DOLFIN_EPS
      totresult = on_boundary and result
  #    if totresult: 
  #      print x, on_boundary, result 
      return totresult
  
  # apply marker1 to positive
  myBoundary = Crowder()
  def doMarking(subdomains,locs,rads,marker):
    for i,x in enumerate(locs): 
      myBoundary.centroid = x
      myBoundary.radius = rads[i]
      myBoundary.mark(subdomains,marker)
  
  
  # mark boudaries 
  mesh = Mesh(xmlName)
  subdomains = MeshFunction("size_t",mesh, mesh.topology().dim()-1)
  subdomains.set_all(0)
  
  #  print case.rads
  for i, case in cases.iteritems():
    alllocs = case.alllocs
    allrads = case.allrads  
    inds= case.inds
    doMarking(subdomains,alllocs[inds],allrads[inds],case.marker)     
  
  #
  ### Solve PB equatio
  #
  #  check that boundaries were marked correctly 
  V = FunctionSpace(mesh,"CG",1)
  m = Function(V)
  # check that markings work 
  for i, case in cases.iteritems():
    bc = DirichletBC(V, Constant(case.V0), subdomains,case.marker)
    bc.apply(m.vector())
  File("markerBoundaries.pvd") << m
  
  
  # apply constants at boundary 
  bcs = []
  for i, case in cases.iteritems():
    bcs.append(DirichletBC(V, Constant(case.V0), subdomains,case.marker))
  #bcs.append(DirichletBC(V, Constant(Vp), subdomains,markerp))
  # Potential: [mV]
  (V,potential ) = pb.PBEngine(mesh,V,subdomains,bcs)
  File("potential.pvd") << potential
  
  scaledPotential = Function(V)
  # PKH - need to check where I am doing the [mV] --> [kT] conversion (e.g 1 kT = 25.6 mV)
  # [mV] *[1/mV] = 1
  scaledPotential.vector()[:] = parms.F_o_RT*potential.vector()[:]
  # [1] *[kT] = [kT] 
  scaledPotential.vector()[:] = parms.kT * scaledPotential.vector()[:]  
  
  print "Potential min/max [kT/e]:"
  ar=np.asarray(scaledPotential.vector()[:])
  print np.min(ar),np.max(ar)
  
  
  # ### Solve smolhomog
  # pass postential to homogenization engine and computer effective diffusivity 
  results= hl.runHomog(fileXML=xmlName,psi=scaledPotential,q=parms.zLig,smolMode=True)
  #results= hl.runHomog(fileXML=xmlName)

  return results 

##
def validation():
  cases = DefineCases_Exp1()
  results = CalcDiffusivityMixedChargeMesh(cases,makeMesh=True,res=10.)

  assert(np.linalg.norm(results.d_eff - np.array([0.79,0.79]))<0.01), "FAILED"
  print "PASS" 

#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

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
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      validation() 
      quit()
 
    elif(arg=="-sel"):
      CalcDiffusivityMixedChargeMesh()
      quit()




  #raise RuntimeError("Arguments not understood")




