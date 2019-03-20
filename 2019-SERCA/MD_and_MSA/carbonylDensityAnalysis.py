import loos

import sys
import math
import PyTraj
import matplotlib.pylab as plt
import numpy as np

class empty:pass

sys.path.append("/home/AD/pmke226/labscripts/parvalbumin/")
import parvfuncs as pf

import binfuncs as bf



## plot binned data
def plotOxyDensity(grid_t,grid_r,zs,title="NONE",ylabel="Mg/O dist [A]"):
  ylimLim = [0,7]
  xlimLim = [0,12]
  nCalCoord = 1.

  scaleByCoordinationNumber=1/nCalCoord
  ar = scaleByCoordinationNumber*np.sum(zs,axis=0)

  # density plot 
  plt.subplot(211)
  plt.title(title)
  plt.xlabel("Time (ns)")
  plt.ylabel("Mg/O dist ($\AA$)")  
  plt.pcolormesh(grid_t/500.,grid_r,ar.transpose())
  plt.colorbar(ticks = (0,1,2,3,4,5,6,7))
  plt.clim([0,7])
  plt.tight_layout()
    
  # histogram plot
  plt.subplot(212)    
  histo=np.sum(ar,axis=1)
  nTs  = np.shape(ar)[1]
  plt.plot(grid_r[0,:],histo/nTs)
  plt.ylim(ylimLim)
  plt.xlim(xlimLim)
  plt.xlabel("Mg/O Distance ($\AA$)")
  plt.ylabel("Carbonyl Density")
  plt.gcf().savefig(title+".png")

  # save data 
  l=np.mean(ar,axis=1) #Same as summing up and averaging over nTs.  See above.
  s=np.std(ar,axis=1)
  j=np.zeros([np.shape(l)[0],3])
  j[:,0]=grid_r[0,:]
  j[:,1]=l
  j[:,2]=s

  np.savetxt(title+".txt",j)
  return ar

## compute density of oxygens relative to Ca position 
def CarbonylDensityAnalysis(case,idx=0,\
                            loopType = "cEF", loopIndices=[88,99], 
                            ionName = "CA", calIdx = 109,
                            binMin=-1,
                            binMax=-1,
                            res = 0.5 # Angstroms
                            ):
    
    ## Make selection
    selAllOs = loos.selectAtoms(case.system,'name=~"O"&&(resid>=%d && resid<=%d)'%(loopIndices[0],loopIndices[1])) # &&segid=="PROA"')
    selCa = loos.selectAtoms(case.system,'name=="%s"&& resid==%d'%(ionName,calIdx))
    #print selAllOs

    ## compute all distancers
    nOs = selAllOs.size()
    nTs = case.trajs[idx].nframes()
    #print "nFrame " , nTs
    allDists = np.zeros((nOs,nTs))
    ts = np.zeros(nTs)
    for i in range(nOs):
        selO = selAllOs.subset(i,1)
        dists = pf.dodist(case,idx,selO,selCa)
        #print np.shape(dists)
        allDists[i,:] = dists[0:nTs,1]
        ts  = dists[0:nTs,0]
      #  print selO

    #myFile = "case0idx%d_dists.txt"%idx
    #np.savetxt(myFile,allDists)
    #myFile = "case0idx%d_ts.txt"%idx
    #np.savetxt(myFile,ts)

    ## get bounds for binning 
    if binMin<0:
        binMin=np.floor(np.min(allDists))
    if binMax<0:    
        binMax=np.ceil(np.max(allDists))
    
    ## do binning for each O-Ca dist
    zs=[]
    for i in range(nOs):
      grid_t,grid_r,z = bf.BinTimeSeries(allDists[i,:],binMin,binMax,res=res,plot=False)#plot=True)
      zs.append(z)
    
    ## display density 
    dummy=plotOxyDensity(grid_t,grid_r,zs,title=loopType)    


def plotOxyDensityAve(ofname1,ofname2,ofname3,title="NONE"):
 o1=np.loadtxt(ofname1,usecols=(0,1,2))
 o2=np.loadtxt(ofname2,usecols=(0,1,2))
 o3=np.loadtxt(ofname3,usecols=(0,1,2))

 data = np.zeros([np.shape(o1)[0],3])
 data[:,0] = o1[:,1]; data[:,1] = o2[:,1]; data[:,2] = o3[:,1]
 mean = np.mean(data,axis=1)
 stdA = np.std(data,axis=1)

 plt.errorbar(x=o1[:,0],y=mean,yerr=std,fmt='g--')
 plt.ylim([-.2,1.4])
 plt.legend(loc=0)
 plt.title(title)
 plt.xlabel('Ca2+-Carbonyl Distance [A]')
 plt.ylabel('Population Density')
 outFile=title+".png"
 plt.gcf().savefig(outFile, dpi=1000)
