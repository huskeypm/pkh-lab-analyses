import numpy as np
import analyzeGotran as anG
import gotranJIT 

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

def GetMonitored(module, ode,tsteps,results,model_params):
    monitored = []
    for expr in sorted(ode.intermediates+ode.state_expressions):
        #print expr.name
        if expr.name not in monitored: # safety
            monitored.append(expr.name)

    monitored_plot = []
    for expr in sorted(ode.intermediates):
        #print expr.name
        if expr.name not in monitored_plot: # safety
            monitored_plot.append(expr.name)

    monitor_inds = np.array([monitored.index(monitor) \
                             for monitor in monitored_plot], dtype=int)
    monitored_get_values = np.zeros(len(monitored), dtype=np.float_)

    # Allocate memory
    plot_values = np.zeros((len(monitored_plot), len(results)))

    for ind, (time, res) in enumerate(zip(tsteps, results)):
        module.monitor(res, time, model_params, monitored_get_values)
        plot_values[:,ind] = monitored_get_values[ monitor_inds ]

    plot_values = np.transpose(plot_values)

    return monitored_plot, plot_values


def runParamsFast(odeName = "shannon_2004.ode",
                  name="out", # if None, returns without writing pickle
                  varDict=None,stateDict=None,
                  dt=0.1,dtn=2000,
                  tsteps = None,  # can explicitly pass in time steps (array)
                  stim_period=False,
                  generateJacobian=False,
                  mxsteps=None,downsampleRate=1,
                  returnDict=dict() # basically a contained for returning results
                  ):
    if stim_period is not False:
        raise RuntimeError("stim period should be set through varDict - not sure why required; fix me")

    params = gotranJIT.init()

    if type( tsteps ) is not np.ndarray:
        params.tstop = dtn
        params.dt = dt
    else:
        print ("Using user-provided time-steps:", np.shape(tsteps))

    if generateJacobian:
        print ("Computing jacobian to save memory")
        params.code.generate_jacobian = True

    # SET Params
    params.parameters = []
    if varDict!=None:
        for key, value in varDict.items():#iteritems():
            params.parameters.extend([key,value])

    # SET states
    if stateDict!=None:
        params.init_conditions=[]
        for key, value in stateDict.items():#iteritems():
            params.init_conditions.extend([key,value])

    ## JIT compile module, simulate
    results,module,tsteps,model_params,ode=gotranJIT.main(odeName, params,tsteps=tsteps)

    # to ensure consistency with old code
    s = results
    t = tsteps
    p = model_params
    s_idx = [state.name for state in ode.full_states]
    p_idx = [param.name for param in ode.parameters]
    assert(len(s_idx)==len(results[0,:])),"Error in state vector. Ask PKH"
    assert(len(p_idx)==len(model_params)),"Error in parameter vector. Ask PKH"

    # get monitored fluxes
    j_idx,j = GetMonitored(module, ode,tsteps,results,model_params)

    if name==None:
        returnDict['data'] = anG.makePackage(p,p_idx,s,s_idx,j,j_idx,t,tsteps)
        return

    if downsampleRate >1:
        sDs,jDs,tDs = downsampleData(s,j,t,downsampleRate)
        red = name
        anG.writePickle(red,p,p_idx,sDs,s_idx,jDs,j_idx,tDs)
    else:
        anG.writePickle(name,p,p_idx,s,s_idx,j,j_idx,t,tsteps)
        returnDict['data'] = anG.makePackage(p,p_idx,s,s_idx,j,j_idx,t)

def downsample(fileName, fileOutName=None,rate=10): 
  dataFull = anG.readPickle(fileName)
   
  p = dataFull['p']
  p_idx = dataFull['p_idx']
  s = dataFull['s']; #sDs = s[::rate,]
  s_idx = dataFull['s_idx']
  j = dataFull['j']; #jDs = j[::rate,]
  j_idx = dataFull['j_idx']
  t = dataFull['t']; #tDs = t[::rate]

  sDs,jDs,tDs = downsampleData(s,j,t,rate)
  
  if fileOutName==None:
    redFileName = fileName.replace(".pickle","_red.pickle")
  else: 
    redFileName = fileOutName 
  anG.writePickle(redFileName,p,p_idx,sDs,s_idx,jDs,j_idx,tDs)

def downsampleData(s,j,t,rate):
  if isinstance(rate,float):
    print ("WARNING: changing rate %f into int")
    rate = np.int( rate ) 
  sDs = s[::rate,]
  jDs = j[::rate,]
  tDs = t[::rate,]
  return sDs, jDs, tDs


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
    namesOnly=False,
    yaml=None):
  # remove pickle
  outBase = outBase.replace(".pickle","")

  # downsample if dt < 1 ms
  if dt < 1.0: # [ms]
    downsampleRateTrial = np.floor(1.0/np.float(dt))
    if downsampleRateTrial > downsampleRate:
      downsampleRate = downsampleRateTrial
      print ("Forcing downsampling rate of ", downsampleRate, " for dt=", dt)


  ### sample param dictionary, can add specific parameters here
  if paramDict==None and yaml is None:
    paramDict = dict()
  elif paramDict is None and yaml is not None:
    paramDict = anG.YamlToParamDict(yaml)
  elif paramDict is not None and yaml is not None:
    raise RuntimeError("choose either, not both")
  print ("{} mM of ATP".format(paramDict["stim_amplitude"]/1e6))

  # stateDict
  #stateDict = None # use defaults for first iter
  # create list of pickle names
  daiters = range(iters)
  pickleNames = [ outBase+"_%d.pickle"%(i+1) for i in daiters ]
  if namesOnly:
    return pickleNames

  # subsequent runs
  for i in range(iters):
      print ("Running iter {} for {} seconds ".format(i+1,dtn) )
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
      mxsteps = 5000
      runParamsFast(odeName=odeName,name=nextName,
                    varDict=paramDict,stateDict=stateDict,dt=dt,dtn=dtn,\
                    mxsteps=mxsteps,downsampleRate=downsampleRate)

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
    print ("t ",tprev)

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
  print ("Downsampling is based on integers/indices; use with caution")
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
    print (pickleNames[0],pickleCatName)
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
      varDict = anG.YamlToParamDict(yaml)
      #varDict = aG.YamlToParamDict(sys.argv[i+1])

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
      print ("PASS!")
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
