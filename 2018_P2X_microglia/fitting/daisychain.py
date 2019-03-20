import numpy as np 
#import analyzeODE as ao
# Grabs stuff from previous run 
#prevOut = "run_G_CaBk1.00_G_NaBk1.00_stim2000_1.pickle"
#prevNum=1 # could probably grab this from pickleName
import analyzeGotran as anG
def InitializeNextInSequence(prevOut,prevNum,downsampleRate):
  # Determine new pickleoutName 
  nextNum=prevNum+1
  nextOut = prevOut.replace("_%d.pickle"%prevNum,"_%d.pickle"%nextNum)

  # Load in prev data 
  data = anG.readPickle(prevOut) 
  si = data['s']
  s_idx = data['s_idx']
  p = data['p']
  p_idx = data['p_idx']


  # create new dict with states/values
  stateDict = dict()
  v = zip(s_idx,si[-1,:]) # grabs  state values from last iteration 
  for i,g in enumerate(v):
    #print "Assigning %f to %s"%(g[1],g[0])
    stateDict[g[0]]=g[1]  

  # create new dict with params/values   
  paramDict = dict()
  v = zip(p_idx,p) # grabs  state values from last iteration 
  for i,g in enumerate(v):
    #print g[0]    
    paramDict[g[0]]=g[1]  
    
  return nextOut,nextNum,stateDict,paramDict   
#s_idx

import runModel as rs
def daisychain(\
    odeName = "shannon_2004.ode",
    dt=0.1,
    dtn=10e3, # elapsed time [ms] per iteration 
    iters = 3,   
    stim_period=1000.,
    mxsteps=None,
    outBase = "run_stim1000",
    stateDict = None,
    paramDict = None,   
    downsampleRate = 1.,
    namesOnly=False):
  # remove pickle
  outBase = outBase.replace(".pickle","")

  # downsample if dt < 1 ms 
  if dt < 1.0: # [ms] 
    downsampleRateTrial = np.floor(1.0/np.float(dt))  
    if downsampleRateTrial > downsampleRate: 
      downsampleRate = downsampleRateTrial
      print "Forcing downsampling rate of ", downsampleRate, " for dt=", dt 
  

  ### sample param dictionary, can add specific parameters here
  if paramDict==None:
    paramDict = dict() 
  else: 
    for key,value in sorted(paramDict.items()):
      if key=="stim_period":
        stim_period = value 
  
  paramDict["stim_period"] = stim_period
  
  # stateDict
  #stateDict = None # use defaults for first iter 
  # create list of pickle names 
  daiters = range(iters) 
  pickleNames = [ outBase+"_%d.pickle"%(i+1) for i in daiters ]
  if namesOnly:
    return pickleNames
  
  # subsequent runs     
  for i in range(iters):
      print "Running iter ", i+1, "for ", dtn, " [ms] at ", stim_period/1e3, "Hz"
      # initialize
      if i==0:
        # first run 
        nextNum = 1
        nextName = outBase+"_%d.pickle"%nextNum
      # second run         
      else:    
        prevName = nextName
        prevNum = nextNum
        nextName,nextNum,stateDict,paramDict = InitializeNextInSequence(\
          prevName,prevNum,downsampleRate)
      
      
      # hack
      #stateDict["V"]=50 works 
      mxsteps = 2000
      rs.runParamsFast(odeName=odeName,name=nextName,
                       varDict=paramDict,stateDict=stateDict,dt=dt,dtn=dtn,\
                       stim_period=stim_period,mxsteps=mxsteps,downsampleRate=downsampleRate)

  ConcatenateTrajs(pickleNames, writeCat=True)

  return pickleNames


# concatenates frajectories 
def concatenateTrajs(pickleNames):
  raise RuntimeError("Antiquated. Use ConcatenateTrajs instead (note different output") 

def ConcatenateTrajs(pickleNames,writeCat=False,downsampleRate=1):
  
  allsi = []
  allji = []
  allt = []
  tprev=0
  for i,pickleName in enumerate(pickleNames):
    data = anG.readPickle(pickleName) 

    si = data['s']
    s_idx = data['s_idx']
    ji = data['j']
    j_idx = data['j_idx']

    p = data['p']
    p_idx = data['p_idx']
    t = data['t']
  
    allsi.append(si)  # probably should pre-allocate 
    allji.append(ji)  # probably should pre-allocate 
    allt.append(t+tprev)    
  
    # update offset  
    tprev += t[-1]  
    print "t ",tprev  

# test 
  #v = [np.arange(5),np.arange(5)+4, np.arange(5)+9]
  #v= np.array(v)
  #print np.ndarray.flatten(v)
  
  ## concatenate times 
  ts = np.array(allt)
  ts = np.array(ts)
  # not general....
  ts = np.ndarray.flatten(ts)
  #print np.shape(ts)
  
  
  ## States 
  # concatenate state values 
  sis = np.array(allsi)
  #print np.shape(sis)
  nSteps = np.prod([np.shape(sis)[0],np.shape(sis)[1]])
  nStates = np.shape(sis)[2]
  #print nSteps
  allsisi = np.zeros([nSteps,nStates])
  # there's a smarter way to hash this out....
  for i in range(nStates):
      sisi = np.ndarray.flatten(sis[:,:,i])    
      allsisi[:,i] = sisi

  ## jis
  # concatenate state values 
  jis = np.array(allji)
  #print np.shape(sis)
  nSteps = np.prod([np.shape(jis)[0],np.shape(jis)[1]])
  nJs     = np.shape(jis)[2]
  #print nSteps
  alljisi = np.zeros([nSteps,nJs])
  # there's a smarter way to hash this out....
  for i in range(nJs):
      jisi = np.ndarray.flatten(jis[:,:,i])    
      alljisi[:,i] = jisi


  #return ts, allsisi, s_idx
  data = dict()
  print "Downsampling is based on integers/indices; use with caution" 
  data['s']    = allsisi[::downsampleRate,]
  data['s_idx']= s_idx 
  data['j']    = alljisi[::downsampleRate,]
  data['j_idx']= j_idx
  data['p']    = p 
  data['p_idx']= p_idx 
  data['t']    = ts[::downsampleRate]

  import re
  if writeCat:
    pickleCatName = re.sub(r'_\d+.pickle', '_cat.pickle', pickleNames[0])
    print pickleNames[0],pickleCatName 
    anG.writePickle(pickleCatName,
                   data['p'],
                   data['p_idx'],
                   data['s'],
                   data['s_idx'],
                   data['j'],
                   data['j_idx'],
                   data['t'])

  return data
  
  
      
  
  

#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################



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
  varDict = dict()
  stateDict = dict()
  odeName = "shannon_2004.ode"
  iters = 3
  dtn = 10e3
  dt = 0.1 # .1 ms default
  stim_period = 1000.
  outBase = "test.pickle"
  downsampleRate=1.

  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if("-var" in arg): # note: can also run yaml 
      # Absolute, not relative values
      varName =sys.argv[i+1]
      varVal =sys.argv[i+2]
      varDict[varName] = np.float(varVal)
   
    if("-yaml" in arg): 
      if bool(varDict): 
        raise RuntimeError("varDict is not empty! Aborting") 
      varDict = aG.YamlToParamDict(sys.argv[i+1])

    if("-state" in arg):
      stateName =sys.argv[i+1]
      stateVal =sys.argv[i+2]
      stateDict[stateName] = np.float(stateVal)
  
    if(arg=="-dSr" or arg=="-downsampleRate" ): 
      downsampleRate= np.float(sys.argv[i+1])

    if(arg=="-dtn" or arg=="-T" ): 
      dtn = np.float(sys.argv[i+1])

    if(arg=="-dt"):
      dt = np.float(sys.argv[i+1])

    if(arg=="-odeName"):
      odeName = sys.argv[i+1]

    if(arg=="-iters"): 
      iters= np.int(sys.argv[i+1])

    if(arg=="-stim_period"): 
      stim_period= np.float(sys.argv[i+1])
  
    if(arg=="-outBase" or arg=="-name"):
      outBase = sys.argv[i+1]

    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      daisychain()
      print "PASS!" 
      quit()

  pickleNames = daisychain(\
    odeName = odeName,
    dt = dt,
    dtn=dtn, # elapsed time [ms]
    iters = iters,
    stim_period = stim_period,
    downsampleRate=downsampleRate,
    outBase = outBase,       
    stateDict = stateDict,
    paramDict = varDict)
    
  





  #raise RuntimeError("Arguments not understood")




