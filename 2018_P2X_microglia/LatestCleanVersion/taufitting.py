import numpy as np 
import matplotlib.pylab as plt 
# nothing 

##
## Fitting
## 

# y = a exp(-t/tau)
def GenData(vars,x):
    amp = vars[0]
    decay = vars[1]

    return amp * np.exp(-x/decay)


# residual (for optimization)
def residual(vars, x, data, eps_data):
    vars[0] = 1 # Fix amplitude, since normalized   
    model = GenData(vars,x)
    return (data-model)/eps_data



from scipy.optimize import leastsq
# Modified to use nonlinear fit where I have zerod out tthe 
# amplitude term 
def FitExp(tsub,caisub,doPlot=False):
  ts = tsub - np.min(tsub)
  ys = caisub - np.min(caisub) + 1e-9
  yn = ys/np.max(ys)

  # linear regression 
  ylog = np.log(yn)
  regression = np.polyfit(ts,ylog, 1)
  aFit = np.exp(regression[1])
  tauFit = -1/regression[0]
  varsFit = [aFit,tauFit]
  print varsFit 

  # nonlinear regression
  eps_data = 1 
  vars = [1,1]      
  out = leastsq(residual, vars, args=(ts, yn, eps_data))
  nlvarsFit = out[0]
  print nlvarsFit

  #print varsFit

  if doPlot:  
    #plt.figure()
    #plt.plot(ts,yn,'b',label="data")
    ##plt.xlim([np.min(ts),np.max(ts)])
    #plt.plot(ts,GenData(varsFit,ts),'g',label="fit")
    #plt.plot(ts,GenData(nlvarsFit,ts),'r',label="nlfit")

    plt.figure()
    plt.plot(ts,caisub,'b',label="data")
    caiest = np.max(ys)*GenData(nlvarsFit,ts) + np.min(caisub)
    plt.plot(ts,caiest,'r',label="nlfit")
    plt.legend()
    #plt.xlim([np.min(ts),np.max(ts)])
  
  return nlvarsFit


##
## Region selection 
## 

# pacing interval 701
# Dafaults to grabbing state data, but can grab fluxes instead
def GetInterval(case,pacingInterval,tstart=8000,  idx=0,getFlux=False,doPlot=False):

    t = case['t']
    if getFlux==False:
      s = case['s']
    else:
      s = case['j']

    if idx!= False: 
      sp = s[1:,idx]
    else: 
      sp = s 
    
    dt = t[1]-t[0]
    span = np.int(pacingInterval/dt)

    si =np.int(tstart/dt)
    sf = si + span


    # assume period starts at max
    sim = np.argmax(sp[si:sf])
    sim+= si
    sfm = sim+span
    sfm = np.argmin(sp[sim:sfm])
    sfm+= sim

    
    if doPlot:
      plt.figure()
      plt.plot(t,sp)
      plt.plot(t[si:sf],sp[si:sf])
      plt.figure()
      plt.plot(t[si:sf],sp[si:sf])
  
      plt.plot(t[sim],sp[sim],"r+")
      plt.plot(t[sfm],sp[sfm],"r+")
      plt.plot(t[sim:sfm],sp[sim:sfm],"k--")
      
    tsub = t[sim:sfm]
    caisub = sp[sim:sfm]
    
    return tsub,caisub


# This doesn't work like I think it does!!
def GetSubRegion(tsteps,cai,subMin,subMax,si):
    #tsteps = t
    # short sims 
    above = np.where(tsteps>subMin)
    #below = np.where(tsteps<2720) # to the base
    below = np.where(tsteps<subMax)
        
    z = np.in1d(above,below)
    inds = above[0][z]

    
    inds2 = inds[np.where(tsteps[inds]>si)]
    tsub = tsteps[inds2]
    #jsub =jsJ_Ca_SL_myo[inds2]/volMyo # divide by volMyo as in .ode code 
    #caisub = s[inds2,idxCai]
    caisub = cai[inds2]
    return tsub,caisub

##
## Tau estimation 
## 
def GetTau(case,pacingInterval,tstart,idxCai,doPlot=False): # subMin=2300,subMax=2600,si=2364):
    ## Find suitable region
    #tsub, caisub = GetSubRegion(tsteps,cai,subMin,subMax,si)
    tsub, caisub = GetInterval(case,pacingInterval,tstart=tstart,idx=idxCai)
    
    #print tsub
    if doPlot:
      plt.figure()
      cais=cai[1:]
      plt.plot(tsteps[subMin:],cais[subMin:],label="test")
      plt.plot(tsub,caisub,"b.", label="sub")
      plt.legend()
    #plt.figure()
    #plt.plot(tsub,caisub,"b.", label="sub")

    ## do fitting 
    fitted = FitExp(tsub,caisub)
    tau = fitted[1]
    #print "tau [ms] ", tau
    return tau

def GetExtreme(tsteps,cai,subMin=2300,subMax=2600,si=2364):
    ## Find suitable region
    tsub, caisub = GetSubRegion(tsteps,cai,subMin,subMax,si)
    
    #print tsub
    doPlot=False
    if doPlot:
      plt.figure()
      cais=cai[1:]
      plt.plot(t[subMin:],cais[subMin:],label="test")
      plt.plot(tsub,caisub,"b.", label="sub")

        
    delCa  =  np.max(caisub)-np.min(caisub)
    print np.min(caisub), np.max(caisub), delCa 
    return delCa  

