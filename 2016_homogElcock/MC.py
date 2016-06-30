#
# Monte Carlo engine for placing non-clasing spheres
#
import matplotlib.pylab as plt
import numpy as np
import miscutil 


# def  generate points/locations that adhere to constraints
def GenLocsRads(
  length = 500., # [A] 
  ptsPerSide = 8, # if False, compute from comDist 
  comDist = False, # if <int>, define pts per side from this 
  ptRad = 10.0,   # [A]
  locScale = 2., # random location scale
  radScale = 1.0,# random radii scale
  eps = 2.,
  keep = 3,
  doPrint=True,
  trials = 1e4,
  inputRads = None, # pass in input set of radii 
  uniform=False): 


  if ptsPerSide==False:
    ptsPerSide = np.floor(length/np.float(comDist))
    print "PtsPerSide %f L/D %f"%(ptsPerSide, length/comDist)



  # determine uniform distributin of points 
  ilocs,rads = SetLocsRads(ptsPerSide,ptRad,length)
  print np.shape(ilocs)
  print np.shape(rads)
 
  # used Passed in radii 
  if inputRads!=None:
    rads = inputRads
    print "Since radii are provided, we are disabling randomization of radii"
    radScale = 0.
  print np.shape(rads)

  # print to file 
  if doPrint:
    temp = np.array(ilocs)
    toPrint = np.zeros((np.shape(temp)[0],3))
    toPrint[:,0:2] = temp
    toPrint[:,2] = rads
    np.savetxt("points.txt",toPrint)

  # for comparison with Comsol results 
  #uniform = True
  if uniform:
    locScale=0;
    radScale = 0.
    keep = 1
    trials = keep 
    
  passingLocs,passingRads = DoMC(ilocs,rads,locScale,radScale,eps,
    chainMode=True,length=length,keep=keep,trials=trials)

  return passingLocs,passingRads


# verify points are within edges of domain ([0,length])
def CheckEdges(nlocs,rads,length,eps=0.01,idx=False):
    rads2 = np.vstack((rads,rads)).T + eps
    test = nlocs*0.

    # check left side
    bottomleft  = nlocs - rads2 
    blClash= (bottomleft < 0)
    test[blClash]=1
    #print bottomleft
    #print blClash

    #right
    topright  = nlocs + rads2 
    trClash=  (topright > length)
    test[trClash]=1.
    #print topright
    #print trClash

    if idx:
      test = np.sum(test,axis=1)
      return np.argwhere(test>0)    
    else:
      return np.any(test>0)                 
   

def SetLocsRads(ptsPerSide,ptRad,length):
    npts = ptsPerSide**2
    rads = np.ones(npts)*ptRad

    # create positions in one D
    pos = np.arange(ptsPerSide)*length/ptsPerSide
    delx = pos[1]-pos[0]
    pos += delx/2.

    # make uniformly distributed in two D
    daones = np.ones(ptsPerSide)
    xs = np.outer(daones,pos)
    ys = xs.T

    # flatten in order to make list of points: locs[N*2]
    xsf = np.ndarray.flatten(xs)
    ysf = np.ndarray.flatten(ys)

    locs = np.vstack((xsf,ysf)).T
    
    return locs,rads

def wrap(a,length):
    b=a.copy()
    
    ind = np.where(a[:,0]<0)
    b[ind,0] = length+a[ind,0]
    ind = np.where(a[:,1]<0)
    b[ind,1] = length+a[ind,1]

    ind = np.where(a[:,0]>length)
    b[ind,0] = a[ind,0]-length
    ind = np.where(a[:,1]>length)
    b[ind,1] = a[ind,1]-length
    
    return b

def MCIter(npts,plocs,rads,locScale,radScale,length,eps):
        # rescale by particle rad
        ls = np.ones(npts)*locScale
        ls/=np.sqrt(rads)
  
        # randomize all locations  
        #rands = np.reshape(np.random.randn(npts,2),(npts,2))
        # randomize one 
        rands = np.zeros((npts,2))
        randNum = np.random.rand(1)*(npts-1)
        randNum = np.int(np.floor(randNum))
        rands[randNum,:] = np.random.randn(2)
        
        rands[:,0]*=ls 
        rands[:,1]*=ls 
  
        # DEBUG rands[:]=locScale  
        tlocs = plocs+rands
  
        # randomize radii   
        rands = radScale*np.random.randn(npts)
        nrads = rads+rands
  
  
        # wrap
        tlocs= wrap(tlocs,length) 
        dx = tlocs - plocs
   
  
        # check   
  #      isClash = miscutil.CheckClash(nlocs,nrads,eps=eps) This failed.
        # Probably should move checkclash here 
        isClash = miscutil.CheckClash(tlocs,nrads,idx=True) 
        #print "asdf",isClash
        isOut = CheckEdges(tlocs,nrads,length,eps=eps,idx=True)
        #print "sdf",isOut
        #isClash = False
        #isOut = False
        #print "Clashing ", isClash
        #print "Out of bounds ", isOut
  
        # reject move 
        #dx[isClash,:]=0
        dx[isOut,:]=0
        isOut=False
        dx[isClash,:]=0
        isClash=False
  
        #print "dx",dx
        # accept rest 
        nlocs=plocs+dx
        #print "nloc",nlocs
        return nlocs,nrads,isClash,isOut 

# Do a bunch of random draws and check for cases that are non-intersecting and within bounds
# chainMode=True, update position based on previous iteration
# length = domain length 
def DoMC(ilocs,rads,locScale,radScale,eps,chainMode=True,keep=5,length=100,trials=1e4):
    npts = np.shape(ilocs)[0]
    
    passingLocs = []
    passingRads = []
    alltrials = []
   
    plocs = ilocs


    ## Equillilibrate 
    for i in range(np.int(trials)):

      nlocs,nrads,isClash,isOut = MCIter(
        npts,plocs,rads,locScale,radScale,length,eps)

      if chainMode:
        plocs = nlocs
  


    ## Get Snapshots 
    #print "pl", plocs
    for i in range(np.int(trials)):
     
      nlocs,nrads,isClash,isOut = MCIter(
        npts,plocs,rads,locScale,radScale,length,eps) 

      if np.any([isClash,isOut]):
            #print "FAIL"
            1
      else:
            #print "SUCCESS"
            passingLocs.append(nlocs)
            passingRads.append(nrads)
      alltrials.append(nlocs)
      
      if chainMode:      
        plocs = nlocs

      #if np.shape(passingRads)[0]>keep:
      #  print i
      #  break
            
    # check        
    if len(passingLocs)==0:
        print "BUMMER: everything failed!! (try different number of points/radii/scale)"
        quit()

    else:
        req=keep 
        keep = np.min([len(passingLocs),keep])
        print "Found %d Keeping %d (of %d requested)" %(len(passingLocs), keep, req)

        
        #passingLocs = passingLocs[0:keep]
        #passingRads = passingRads[0:keep]

        #passingLocs = passingLocs[-keep:]
        #assingRads = passingRads[-keep:]

        #arr = np.random.shuffle(np.arange(len(passingLocs)))

        # grab a few at random 
        arr = np.arange(len(passingLocs))
        np.random.shuffle(arr)
        inds = arr[0:keep]
        #print arr
        print inds
        l=[]
        r=[]
        for i in inds:
          l.append(passingLocs[i])    
          r.append(passingRads[i])    
        passingLocs=l
        passingRads=r
      
  
    return passingLocs,passingRads


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

def Test2():
  length = 600
  radGyr = 30.
  minEdgeDist = 10.
  locScale=2.0  
  radScale=2.0 
  keep =4
  comAvg = 2*radGyr+minEdgeDist

  passingLocs, passingRads = GenLocsRads(
     length=length,
     keep=keep,
     ptsPerSide=False,
     ptRad = radGyr,
     comDist = comAvg,
     locScale = locScale, # these needed to be determined from the protein data
     radScale = radScale,# these needed to be determined from the protein data.
     trials=1e6
      )

  for i,nlocs in enumerate(passingLocs):
  #for i,nlocs in enumerate(alltrials):
    nrads = passingRads[i]
    plt.figure()
    miscutil.circles(nlocs[:,0],nlocs[:,1],nrads)

    plt.xlim([0,length])
    plt.ylim([0,length])
    plt.axes().set_aspect('equal')

    plt.gcf().savefig("test%d.png"%i,dpi=300)


def doit():
  length = 100.  # [A] 
  ptsPerSide = 3
  ptRad = 4.0    # [A]
  locScale = 2.  # random location scale
  radScale = 1.0 # random radii scale
  eps = 1.
  
  locs,rads = SetLocsRads(ptsPerSide,ptRad,length)
  passingLocs,passingRads = DoMC(locs,rads,locScale,radScale,eps,keep=3,length=length)
  
  for i,nlocs in enumerate(passingLocs):
  #for i,nlocs in enumerate(alltrials):
    nrads = passingRads[i]
    plt.figure()  
    miscutil.circles(nlocs[:,0],nlocs[:,1],nrads)
  
    plt.xlim([0,length])
    plt.ylim([0,length])
    plt.axes().set_aspect('equal')

    plt.gcf().savefig("test%d.png"%i,dpi=300)





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
    if(arg=="-test2"):
      Test2()
      quit()
  





  raise RuntimeError("Arguments not understood")




