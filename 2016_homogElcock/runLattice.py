import numpy as np 
import matplotlib.pylab as plt
from scipy  import meshgrid

import cPickle as pickle

class empty:pass

import sys
root = "/home/AD/pmke226/sources/"
sys.path.append(root+"/smolhomog/example/validation/")
import buildMesh as bm
import miscutil as mu

sys.path.append(root+"/smolfin/")
sys.path.append(root+"/homogenization/")
import homoglight as hl

def HSBoundCyl(phi):
    return phi*2/(3-phi)

# replaces newrandomizer-150516
def RunMCMesh():
  # List of cell sizes 
  #unit_cell_size=np.array([300,400,500,600,700,800]) # [A] 
  unit_cell_size=np.array([300])
  nCs=np.shape(unit_cell_size)[0]



#import runLattice as rL
# For testing the effect of lattice defects on macroscale properties 
def RunDefectsLattice(storeResults=True,nDefects=-1):
  nPerSide = 5
  lcell = 200 # 
  radGyrAvg=30. 
  instances=5
  if nDefects<0:
    defects = np.array([0,1,2,3,4,5]) # randomized 
    nDefectTests=np.shape(defects)[0]
  else: 
    instances=1
    defects = np.array([nDefects])
    nDefectTests = 1
  
  defectis=np.zeros([nDefectTests,instances],dtype=int)
  dxis=np.zeros([nDefectTests,instances])
  for i in np.arange(nDefectTests):
      for j in np.arange(instances):
        locs = makeRegLattice(lcell,nPerSide, (nPerSide*nPerSide-defects[i]))
        # debg dxi = i+np.random.rand(1)  
        dxi,dyi,resu = HomogLattice(locs=locs,radGyrAvg=radGyrAvg)
        defectis[i,j] = defects[i]      
        dxis[i,j]=dxi

   # store 
  if storeResults:
    data1 = {'defects':defects,'defectis':defectis,'dxis':dxis}
    output = open('defect.pickle', 'wb')
    pickle.dump(data1, output)
    output.close()

        

# Generate 2D mesh from list of 2D points 
def HomogLattice(locs,radGyrAvg=1.,\
                 res = 3.0,   # Tet resolution [Same units as locs]
                 useRandom=False,\
                 computeMargin=True,   
		 boxMin=[1,1],boxMax=[0,0],
		 rads=None, # list of radii 
                 gmshName="out.geo",   
                 pb=None,   # do poisson boltzmann solve ,
                 skipDolfin=False,
                 skipMeshGen=False):
    ## special options 
    if useRandom:
      raise RuntimeError("No longer supported") 
    #    # might have to run a couple of times if intersecting
    #    print "Generating random points" 
    #    pts = 10
    #    locs = 100*np.random.rand(pts,2)
    #    rads = 2+np.random.rand(pts)
    #    marg = 1.  

    pts =  np.shape(locs)[0]
    #print "locs: ",locs
    #print "pts:", pts
    if rads==None:    
        rads = np.ones(pts) * radGyrAvg     

    if computeMargin:
        #rads = 2+np.random.rand(np.shape(locs)[0])
        marg = 0.5*radGyrAvg
    
        maxRad = np.max(rads)
        boxMin = np.min(locs,axis=0) - marg - maxRad
        boxMax = np.max(locs,axis=0) + marg + maxRad
  
        
    ## Here's the meat 
    # Gen mesh 
    print "Creating mesh with %d points " % pts
    if skipMeshGen:
      xmlName = "out.xml"
    else:
      xmlName = bm.makeLattice(locs,rads,boxMin=boxMin,boxMax=boxMax,\
                               res=res,fileName=gmshName)


    if skipDolfin:
      return [-1,-2,-3]

    # add in pb component 
    if pb!=None:
      parms = pb.parms
      # assume boundary potential can be determined via Grahame eqn 
     # if parms.V0!= None and parms.sigma!=None:
     #   raise RuntimeError("You cannot define V0/sigma simultaneously") 
      # [mV]
      if parms.V0 == None:
        parms.V0 = pb.Grahame(parms.sigma,parms.ionC)

      # Solve PB eq subject to boundary bc
      from dolfin import Mesh,Function 
      mesh = Mesh(xmlName)
      (V,potential)= pb.SolvePoissonBoltzmann(
          mesh,boundaryPotential=parms.V0,boundaryCondition="crowders",
          outerZero=False, # do not set potential to zero at outer boundary 
          outName=xmlName.replace(".xml","_pb.pvd"))
      mu.Plot2DSolution(mesh,potential,fileName=xmlName.replace(".xml","_pb.png"),dpi=100,clim=[parms.V0,0]) 

      # convert to PMF 
      scaledPotential = Function(V)
      # PKH - need to check where I am doing the [mV] --> [kT] conversion (e.g 1 kT = 25.6 mV)
      scaledPotential.vector()[:] = parms.F_o_RT*potential.vector()[:]
      scaledPotential.vector()[:] = parms.kT * scaledPotential.vector()[:]  # (z=+,psi=+) --> [mV]

      results= hl.runHomog(fileXML=xmlName,psi=scaledPotential,q=parms.zLig,smolMode=True)
    # no PB (neutral 
    else:
      results = hl.runHomog(fileXML=xmlName)
    
    # Get analytical bounds from HS theory
    dHS = HSBoundCyl(results.volFrac)
    results.dHS = dHS
    print "Deff/HS", results.d_eff, dHS
    
    return results.d_eff[0],results.d_eff[1],results

# selct points within bounds 
def makeIrregLattice(locs,boundsMin,boundsMax):
  idx = np.nonzero( (locs[:,0] > boundsMin[0]) & (locs[:,0] < boundsMax[0]) &\
                      (locs[:,1] > boundsMin[1]) & (locs[:,1] < boundsMax[1]) )
  #print locs[idx,]
  locsBounded =(locs[idx,])[0]
  return locsBounded

# #### Define unit cell
# 
# $l_{cell} = \sqrt{\frac{\pi R_g^2}{\rho_a}}$
def makeRegLatticeFromIrreg(locsBounded,boundsMin,boundsMax,radGyrAvg):

  # compute area density based on estimate of average radius of gyration 
  areaDomain = np.prod(boundsMax-boundsMin)
  nParts = np.shape(locsBounded)[0]
  if nParts<1:
    return []
  areaParts = nParts * np.pi * radGyrAvg**2
  rhoa =(areaParts/areaDomain)
  print "Area density %f " % rhoa



  # <codecell>
  lcell = np.sqrt((np.pi*radGyrAvg**2)/rhoa)
  print lcell

  # overestimate number of regular points and make regular lattice 
  ceil = np.ceil(np.sqrt(nParts))
  nPerSide = np.int(ceil)
  locsReg= makeRegLattice(lcell,nPerSide,nParts=nParts)
  return locsReg  
    
def makeRegLattice(lcell,nPerSide,nParts=-1):    
  x,y = meshgrid(np.arange(nPerSide),np.arange(nPerSide))
  #latticeLocs = np.outer(np.arange(nPerSide),np.arange(nPerSide))
  locsReg = np.asarray( zip( np.ndarray.flatten(x), np.ndarray.flatten(y)) )
#  print "locsReg: ", locsReg
  locsReg*=lcell
#  locsReg*=100
#  print "locsReg: ", locsReg
  locsReg+= np.ones(2) * lcell/2.
  #print "locsReg: ", locsReg

  # select only subset
  if nParts>0:
    import random
    idxs =  random.sample(np.arange(nPerSide**2),nParts)
    locsReg = locsReg[idxs,:]
  
  #print np.shape(locsReg)
  return locsReg

# computes effective diffusion parameters for a rnage of unit cell sizes, based on MC generated crowder positions
def RunLatticesMC(radGyrAvg,comAvg,unit_cell_sizes, \
                  randLocs=True,randGyrs=True,debug=False,fileName=None):
  import MC
  locScale = 0.
  radScale = 0.
  eps = 2.
  keep=3
     
  if randLocs:   
    locScale = 2.  # random location scale
  if randGyrs:        
    radScale = 5.0 # random radii scale
             
 
     
  nCs = np.shape(unit_cell_sizes)[0]
  xDeffs=np.zeros((nCs,keep))
  yDeffs=np.zeros((nCs,keep))
 
  dCLave=np.zeros(nCs)
  dCLstd=np.zeros(nCs)
  dCLvolFracAve=np.zeros(nCs)    
 
  randScale = 1.0 # (parameter controlling size of distribution)
  radGyrStdv = 1.
  
  for i in range(len(unit_cell_sizes)):
     length = unit_cell_sizes[i]   # dimension of a 2D cell (in Angstroms [A])         
     
     passingLocs, passingRads = MC.GenLocsRads(
       length=length, ptsPerSide=False, # use COM avg instead
       ptRad=radGyrAvg,comDist=comAvg,
       locScale=locScale,radScale=radScale,keep=keep)
 
     # This should be done after the MC gen. 
     #    rads_tot=np.append(rads_tot,nrads)
     #    plt.figure()
     #    hist,bins,ignored=plt.hist(nrads,normed=True)
     #    plt.plot(bins,mlab.normpdf(bins,scale,10),linewidth=2,color="r")
     
 #if 0:            
     #rads_ar=np.array(rads_tot)
     #plt.figure()
     #hist,bins,ignored=plt.hist(rads_ar,normed=True) #get hist
     #plt.plot(bins,mlab.normpdf(bins,scale,10),linewidth=2,color="r")
  # plt.plot(bins,1/(10*np.sqrt(2*np.pi))*np.exp(-(bins-scale)**2/(2*10**2)),linewidth=2,color="k")
  # plt.gcf().savefig("/net/share/cesc235/homog/gauss_dist.png")
     #(mus,sigma)=norm.fit(rads_ar) 
     #print "Stats"
     #print mus, sigma
 
     
     if debug:
         continue
         
     # do homog         
     dxs,volFracs = runLatticeNonUnif(passingLocs,passingRads,length,k=i)        
     dCLave[i] = np.mean(dxs)
     dCLstd[i] = np.std(dxs)
     dCLvolFracAve[i] = np.mean(volFracs)
  
  # store for later retrieval 
  if fileName!=None:   
    outar = np.zeros([np.shape(dCLave)[0],4])
    outar[:,0] = unit_cell_sizes
    outar[:,1] = dCLave
    outar[:,2] = dCLstd
    outar[:,3] = dCLvolFracAve
    #print outar
    np.savetxt(fileName,outar)  
     
  return dCLave, dCLstd, dCLvolFracAve    
    


# runs homog on a lattice with non-uniformly distributed crowder locations
# and radii 
def runLatticeNonUnif(
  passingLocs,# list of arrays of locations
  passingRads,# list of arrays of radii
  length, # edge length  
  k=0, # case label
  doPlot = False,
  skipMeshGen = False,
  debugMesh = False, # for print geo files only 
  res=10.,
  pb=None     # if defined, we do the Poisson Boltzmann stuff 
  ):

  dxs = []
  volFracs = []

  for i,nlocs in enumerate(passingLocs):
  #for i,nlocs in enumerate(alltrials):
    nrads = passingRads[i]
    if doPlot:   
        #plt.figure()  
        # not sure why this isn't working now  
        #miscutil.circles(nlocs[:,0],nlocs[:,1],nrads)aa
        #print nlocs
        #print nrads  
  
        #plt.xlim([0,length])
        #plt.ylim([0,length])
        #plt.axes().set_aspect('equal')
        continue
  
          
    boxMin = [0,0]
    boxMax = [length,length] 
          
    pts =  np.shape(passingLocs)[0]
    dxr,dyr,results = HomogLattice(nlocs,\
      res = res,      
      rads=nrads,      
      computeMargin=False,
      boxMin=boxMin,boxMax=boxMax,
      gmshName = "temp_Draw%d_Cell%d.geo"%(i,k),
      skipMeshGen = skipMeshGen,
      skipDolfin=debugMesh,
      pb=pb                
    )
    if debugMesh==False:
      dxs.append(dxr)
      volFracs.append(results.volFrac)  
      print "volume Fraction ", volFracs

  return dxs,volFracs


# Runs homog on a 2D lattice of points 
# given a list of 2D points (locs), this routine
# selects all points within boundsMax-boundsMin
# Should be renamed
def runLattice(locs,boundsMin,boundsMax,radGyrAvg=1.,doPlot=False):
  # name case 
  name = "cell_%.3d_%.3d_%.3d_%.3d" % \
         (boundsMin[0],boundsMin[1],boundsMax[0],boundsMax[1])
  

  locsIrreg= makeIrregLattice(locs,boundsMin,boundsMax)

  # <markdowncell>
  locsReg =  makeRegLatticeFromIrreg(locsIrreg,boundsMin,boundsMax,radGyrAvg)

  if doPlot: 
    plt.figure()
    plt.scatter(locsReg[:,0],locsReg[:,1],edgecolor='red',facecolor="none")
    plt.scatter(locsIrreg[:,0],locsIrreg[:,1])
    plt.gcf().savefig(name+".png")
  
  # ## Run homogenization on different lattice arrangements 
  dxr,dyr,resultsr= HomogLattice(locs=locsReg,radGyrAvg=radGyrAvg)
  dxi,dyi,resultsi= HomogLattice(locs=locsIrreg,radGyrAvg=radGyrAvg)          
  
  # <codecell>
  resultsr.dhs = 0.071
  resultsi.dhs = 0.071
  return [dxr,dyr,resultsr.dhs],[dxi,dyi,resultsi.dhs]
  
  
#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#
def doit():
    locs = np.array([[50,50],[100,100]])
    boundsMin = np.array([0,0])
    boundsMax = np.array([370,370])
    runLattice(locs,boundsMin,boundsMax,radGyrAvg=30.,plot=False)


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
      doit()      
      exit()
    if(arg=="-defects"):
      RunDefectsLattice()
      exit()
    if(arg=="-defectsTest"):
      RunDefectsLattice(storeResults=False,nDefects=1)
      exit()
  





  raise RuntimeError("Arguments not understood")




