import numpy as np
import matplotlib.pylab as plt
import loos
import PyTraj



import parvfuncs as pf



cols = ["-","-.","--"]

def seldist(dists,fr,sel1,sel2):
  c1 = sel1.centroid()
  c2 = sel2.centroid()
  diff = c1-c2
  distance = diff.length()
  dists[fr,1] = distance    

def dodist(case,run,sel1,sel2):
    dists = np.zeros((case.trajs[run].nframes(),2))
    dists[:,0] = np.arange(case.trajs[run].nframes())*case.frame_to_time
    
    # instantiate a new iterator 
    ptraj = PyTraj.PyTraj(case.trajs[run], case.system)
    
    fr=0
    for frame in ptraj:
      seldist(dists,fr,sel1,sel2)  
      #seldist(distschain2,fr,sel1chain2,sel2chain2)  
      fr+=1
    print fr        
    return dists

# barplot
def BarPlotMeans(cases,runs,title="",root="./"):
    plt.figure()
    fig,ax = plt.subplots()
    width = 0.25
    for j, run in enumerate(runs):
        means = []
        stds = []
        nsamples = []
        idxs = np.arange(len(cases))
        names = []
        barcolors = []
        for i, case in enumerate(cases):
          means.append(case.means[j])
          stds.append(case.stds[j])  
          nsamples.append(case.nsamples[j])    
          names.append(case.name)            
          barcolors.append(case.col)
            
            
        # std error
        stds = np.array(stds)
        nsamples = np.array(nsamples)
        #conf = 1.96 * stds/np.sqrt(nsamples)
        #barlist = ax.bar(idxs+j*width,means,yerr=conf,width=width)
        #print stds
        barlist = ax.bar(idxs+j*width,means,yerr=stds,width=width,ecolor='k')
        for i, case in enumerate(cases):
          barlist[i].set_color(barcolors[i])

    plt.title(title)
    ax.set_xticks(idxs+1.5*width)      
    ax.set_xticklabels( names )
    plt.gcf().savefig(root+"barplot"+title+".png")

def PlotCalphaDist(cases,runs,res1,res2,\
                   adjylim=False, ylim=[0,0],\
                   title="",root="/net/share/anku223/PV/"):
    ## compute dists
    ## plot trajectory 
    caseName = "unk"
    plt.figure()
    for i,case in enumerate(cases):
        sel1 = loos.selectAtoms(case.system,'name == "CA" && resid == %d' % res1)  
        sel2 = loos.selectAtoms(case.system,'name == "CA" && resid == %d' % res2)          

        case.means = np.zeros(np.int(len(runs)))
        case.stds = np.zeros(np.int(len(runs)))
        case.nsamples = np.zeros(np.int(len(runs)))
        for j, run in enumerate(runs):
          idx = np.int(run)-1 # since not zero-indexes    
          dists = dodist(case,idx,sel1,sel2)
          plt.plot(dists[:,0],dists[:,1],case.col+cols[j],label="%s/%s"%(case.name,run))   

          metric = dists[:,1]
          case.means[j] = np.mean(metric)
          case.stds[j] = np.std(metric)
          case.nsamples[j] = np.shape(metric)[0]
          #print "mean %f/std %f" % (case.means[j],case.std[j])   

    #case.
    plt.xlabel("time [ns]")
    #plt.legend(bbox_to_anchor=(1.45, 1), borderaxespad=0.)
    plt.legend(loc=0, ncol=2)
    if adjylim:
      plt.ylim(ylim)
    plt.title(title)

    plt.gcf().savefig(root+caseName+"_res%dres%d_CaDist.png" % (res1,res2))

    ## compute bar plot 
    BarPlotMeans(cases,runs,title=title,root=root) # "CalphaDist")



# Helix posn/tilt wrt starting config

def helixInfo(case,h,run=0,skip=0,stride=1):

    calphas = loos.selectAtoms(case.system, 'name == "CA"')
    #calphas = loos.selectAtoms(case.system, selectionString)

    
    # Align based on calphas, but iterate over all atoms 
    atraj = PyTraj.PyAlignedTraj(case.trajs[run],case.system,\
                                 alignwith = calphas, skip = skip, stride=stride)
    selectionString = 'name == "CA" && resid >= %d && resid <= %d' % (h[0],h[1])
    sel1 = loos.selectAtoms(case.system,selectionString)

    # data vectors 
    totFrames = case.trajs[run].nframes()
    nFrames = np.ceil((totFrames-skip)/np.float(stride))
    frameRate = case.frame_to_time* stride # WILL APPEAR INCORRECT IF AUTOLOAD USED FOR MODULE  
    angles = np.zeros((nFrames,2))
    #print np.shape(angles)
    #print np.shape(np.arange(nFrames))
    angles[:,0] = np.arange(nFrames)*frameRate
    centroids = np.zeros((nFrames,3))
    centroidDists = np.copy(angles)

    # iterate 
    #DBGv0 = np.array([1,0,0])
    fr=0
    for frame in atraj:    

      # angles 
      pca = sel1.principalAxes()
      vi = pca[0]
      if vi.z() < 0:
        vi *= -1 
      vi = np.array(vi)
      nvi = np.linalg.norm(vi)

      # store 
      if fr==0:
        c0 = sel1.centroid()
        v0=vi
        nv0=nvi
        
      #DBGvi = np.array([1,0,0])
      cosine = np.dot(v0,vi) / (nv0*nvi)
      angles[fr,1]= cosine

      # centroids
      #print np.shape(centroids[fr,:])
      c = sel1.centroid()
      centroids[fr,:]=[c[0],c[1],c[2]]        
       
      # get distance 
      centroidDists[fr,1]=np.linalg.norm(c0-sel1.centroid())
      #print c0, sel1.centroid(),centroidDists[fr,1]
      #print v0, vi

      fr+=1

    # compute angles 
    cosines = np.round(angles[:,1],decimals=5)    
    angles[:,1] = np.arccos(cosines)*180./np.pi    
    
    return centroids,centroidDists,angles


def doHelixPlot(cases,runs,h,helixName,stride=1,skip=0):
    print "HERE"
    ## Get helix information 
    # print "stride/skip", stride, " ", skip
    for i,case in enumerate(cases):
        case.centroids = [None] * len(runs)
        case.centroidDists = [None] * len(runs)
        case.angles = [None] * len(runs)
        case.means   = [None] * len(runs)
        case.stds    = [None] * len(runs)
        case.nsamples= [None] * len(runs)
        for j, run in enumerate(runs):
          centroids, centroidDists, angles = helixInfo(case,h,run=j,stride=stride,skip=skip)#, stride=10)#,skip=7000)
          case.centroids[j]  = np.copy(centroids)        
          case.centroidDists[j]  = np.copy(centroidDists)        
          case.angles[j] = np.copy(angles)

          #metric = angles[:,1]       
          #title ="Angles"
          metric = centroidDists[:,1]       
          title ="Centroid Dists"
          case.means[j] = np.mean(metric)
          case.stds[j] = np.std(metric)
          case.nsamples[j] = np.shape(metric)


    ## Plot helix information 
    ax=plt.subplot(2,1,1)
    plt.title(helixName + " (%d-%d)"% (h[0],h[1])+ " Angle")
    for i,case in enumerate(cases):
        for j, run in enumerate(runs):
            angles = case.angles[j]
            #print case.col, " ", cols[j]
            #plt.plot(angles[:,0],angles[:,1],label="%s/%s"%(case.name,run)) 
            plt.plot(angles[:,0],angles[:,1],case.col+cols[j],label="%s/%s"%(case.name,run))   
    ax.set_xticklabels([])
    #plt.xlabel("Time [ns]")
    plt.ylabel("Angle [deg] wrt t=0")
    plt.legend(bbox_to_anchor=[1.4,1])

    plt.subplot(2,1,2)
    plt.title(helixName + " CentroidDist wrt t=0")
    for i,case in enumerate(cases):
        for j, run in enumerate(runs):
            centroidDists = case.centroidDists[j]
            #plt.plot(centroidDists[:,0],centroidDists[:,1],label="%s/%s"%(case.name,run))   
            plt.plot(centroidDists[:,0],centroidDists[:,1],case.col+cols[j],label="%s/%s"%(case.name,run))   

    plt.xlabel("Time [ns]")
    plt.ylabel("Distance [Ang]")
    plt.gcf().savefig("helixPlots_"+helixName+".png")
    
    BarPlotMeans(cases,runs,title=title)

def CalcInterhelicalDist(cases,runs,h1,h2,skip=0,stride=1):
    # get distances based on prior calcs 
    for i,case in enumerate(cases):
        case.interHelicalDists = [None] * len(runs)
        case.times = [None] * len(runs)
        
        case.means = np.zeros(np.int(len(runs)))
        case.stds = np.zeros(np.int(len(runs)))
        case.nsamples = np.zeros(np.int(len(runs)))        
        for j, run in enumerate(runs):
            #case.times[j] = np.arange(case.trajs[j].nframes())*FRAME_TO_TIME
            totFrames = case.trajs[np.int(run)-1].nframes()
            nFrames = np.ceil((totFrames-skip)/np.float(stride))
            frameRate = case.frame_to_time* stride           
            case.times[j] = np.arange(nFrames)*frameRate
            
            # sel 1 
            if h1=="HA":
              sel1 = case.hAcentroids[j]
            elif h1=="HC":
              sel1 = case.hBcentroids[j]
            elif h1=="HE":
              sel1 = case.hEcentroids[j]
            else:
              raise RuntimeError("Need to add")
                
            # sel2                 
            if h2=="HD":
              sel2 = case.hDcentroids[j]
            elif h2=="HF":
              sel2 = case.hFcentroids[j]
            else:
              raise RuntimeError("Need to add")
                
            case.interHelicalDists[j] = np.linalg.norm(sel1-sel2,axis=1)
            metric = case.interHelicalDists[j]
            case.means[j] = np.mean(metric)
            case.stds[j] = np.std(metric)
            #print np.shape(metric)
            case.nsamples[j] = np.shape(metric)[0]
            
    
    # Now plot             
    name = h1+h2
    plt.figure()
    plt.title("%s/%s dist" % (h1,h2)) 
    for i,case in enumerate(cases):
        case.means = [None]*len(runs)
        for j, run in enumerate(runs):
            centroidDists = case.interHelicalDists[j]
            case.means[j] = np.mean(centroidDists)
            plt.plot(case.times[j],centroidDists,case.col+cols[j],label="%s/%s"%(case.name,run))   
    plt.xlabel("Time [ns]")
    plt.ylabel("Distance [Ang]")
    plt.gcf().savefig("helixPlots_"+name+".png")
    plt.legend(bbox_to_anchor=(1.3,1.0))
    plt.gcf().savefig("centroidDist_"+name+".png")
            
    BarPlotMeans(cases,runs,title="Interhelical dists")



def ResHelixDist(cases,runs,resid,hIdx):
  plt.figure()
  for i,case in enumerate(cases):
    sel1 = loos.selectAtoms(case.system,'name == "CA" && resid == %d' % resid)
    sel2 = loos.selectAtoms(case.system,'name == "CA" && resid >= %d && resid <= %d' % (hIdx[0],hIdx[1])) 
  
    case.means = np.zeros(np.int(len(runs)))
    case.stds = np.zeros(np.int(len(runs)))
    case.nsamples = np.zeros(np.int(len(runs)))
      
    for j, run in enumerate(runs):
      idx = np.int(run)-1 # since not zero-indexes    
      dists = pf.dodist(case,idx,sel1,sel2)
      
      metric = dists[:,1]
      case.means[j] = np.mean(metric)
      case.stds[j] = np.std(metric)
      case.nsamples[j] = np.shape(metric)[0]
      
      plt.plot(dists[:,0],dists[:,1],case.col+cols[j],label="%s/%s"%(case.name,run) )

  title="Res%d/helix (%d-%d) distance" % (resid,hIdx[0],hIdx[1])
  plt.title(title)   
  plt.legend(bbox_to_anchor=[1.5,1])

  plt.figure()
  BarPlotMeans(cases,runs,title=title)

def ElectroEnergy(fname="NONE"):
#  fname="/home/AD/cesc235/labscripts/efhands/WTholo_pEF_chainA_run0.txt"
  dist,height,stdev=np.loadtxt(fname,usecols=(0,1,2),unpack=True)
#print dist
#print height
#print stdev

  qCa=2
  qCO=-0.70
  
  num=np.shape(dist)[0]
  prob=np.zeros(num)
  probability=np.array(prob)
  h=np.array(height)
  d=np.array(dist)

  Nold=0
  Nsum=0

  i=0
  for i in range(num):
    probability[i]=6*h[i]
    if d[i] != 0:
        N=(qCO*probability[i])/d[i]
        Nsum=Nold+N
#        print i, d[i], N
#        print Nsum
        Nold=Nsum
  Eele=Nsum*qCa
#  print Nsum
  return Eele 
 # return Nsum

def plotAveOxyDensity(ofname1,ofname2,ofname3,title="NONE"):
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
 plt.gcf().savefig(outFile)
