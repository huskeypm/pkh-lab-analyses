"""
For performing parameter sweeps and running model with dictionaries of parameter values 
"""

# Revisions 
# Changed default stimulation to 701
import numpy as np
import runner 
import downSamplePickles
runner.init()
#idxNCX = runner.model.monitor_indices("i_NaCa")
#import analyzeODE as ao
import analyzeGotran as anG
import gotranJIT

class empty:
    def __init__(self,p,s,ts,js):
      self.p = p
      self.s = s
      self.ts = ts
      self.js = js

###
###
###
# currently provides a mechanism for creating a series of 
# time steps with differing resolutions
def TimeSteps(
  dtFine  =0.1,
  dtnFine = 3e3, # 0..dtnFine with dt=dtFine
  dtCoarse=10.,  # dtnFine+dtFine..dtnCoarse with dt=dtCoarse
  dtnCoarse=9e3
  ): 

  ts=0
  tstepsFine =   np.linspace(ts, dtnFine, (dtnFine-ts)/dtFine+1)
  ts=dtnFine+dtFine
  tstepsCoarse = np.linspace(ts, dtnCoarse, (dtnCoarse-ts)/dtCoarse+1)
  #print np.shape(tstepsFine),tstepsFine[0],tstepsFine[-1]
  #print np.shape(tstepsCoarse),tstepsCoarse[0],tstepsCoarse[-1]
  tsteps = np.concatenate([tstepsFine,tstepsCoarse])
  return tsteps 


###
### JIT support  
###

# Shamelessly snagged from gotranrun.py 
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

#tstop=20000 # works 
#tstop=180000 # kicks the bucket after 130 s
#dt = 0.1
#stim_period = 1000.  # works (after adjusting odeint params)
#stim_period = 500.
# WARNING: here that parameters are SET, not rescaled, in contrast to runParams function
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
    print "Using user-provided time-steps:", np.shape(tsteps) 

  if generateJacobian:
    print "Computing jacobian to save memory"
    params.code.generate_jacobian = True
  #params.downsampleRate = downsampleRate

  # SET Params
  #params.parameters = ['stim_period',stim_period]  
  params.parameters = []  
  if varDict!=None:
    for key, value in varDict.iteritems():
      params.parameters.extend([key,value])  
  #print params.parameters
  # SET states 
  if stateDict!=None:   
    params.init_conditions=[]
    for key, value in stateDict.iteritems():
      #print key, value 
      params.init_conditions.extend([key,value])  
    
  ## JIT compile module, simulate  
  #ks=25. # default  [1/ms]
  #KSRleak=5.348e-6 # default [1/ms]
  #params.parameters = ['ks',ks,'KSRleak',KSRleak]  
  #print params.parameters    
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
  #print "s", "j", "t", np.shape(s), np.shape(j), np.shape(t)
  #print "s_for_real", s
  #print "ji", len(j_idx)   
  
  #print "time", np.min(t), np.max(t)
  if name==None:
    returnDict['data'] = anG.makePackage(p,p_idx,s,s_idx,j,j_idx,t,tsteps)
    return 

  if downsampleRate >1:     
      sDs,jDs,tDs = downSamplePickles.downsampleData(s,j,t,downsampleRate)
#      print "jds", np.shape(jDs)
#     print "name ", name 
      red = name 
      anG.writePickle(red,p,p_idx,sDs,s_idx,jDs,j_idx,tDs)
  else:
      anG.writePickle(name,p,p_idx,s,s_idx,j,j_idx,t,tsteps)
      returnDict['data'] = anG.makePackage(p,p_idx,s,s_idx,j,j_idx,t)

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
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg += "  %s -sweep nameVar startVal endVal incrVal " % (scriptName)
  msg += " \n or\n "
  msg += "  %s -var nameVar val  " % (scriptName)
  msg += """
  
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
  pi = runner.model.p
  stim = 1000 # [ms] 
  dt = 1. # [ms] 
  name="out"
  deltaT = 10000 # [ms] 
  downsampleRate = 1  
  sweep = False
  useJIT=True # There shouldn't be a compelling reason to set this to false 
  varDict = dict()              
  odeName = "shannon_2004.ode"
  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if("-odeName" in arg): 
      odeName = sys.argv[i+1]
      # wont work for useJIT=False
    if("-var" in arg):
      varName =sys.argv[i+1] 
      varVal =sys.argv[i+2] 
      varDict[varName] = np.float(varVal)
    if("-sweep" in arg):
      varName =sys.argv[i+1] 
      varVals =sys.argv[(i+2):(i+5)] 
      varDict[varName] = varVals
      sweep=True
      
    if(arg=="-dt"):
      dt=np.float(sys.argv[i+1])
      
    if(arg=="-T"):
      deltaT=np.float(sys.argv[i+1])
    #if(arg=="-downsampleRate"):
    if(arg=="-dSr" or arg=="-downsampleRate" ):
      downsampleRate = np.float(sys.argv[i+1])
    if(arg=="-name"):
      name=sys.argv[i+1] 
    if(arg=="-stim"):
      stim=sys.argv[i+1] 

    if(arg=="-jit"):
      useJIT = True

  # execute
  if sweep:
    if useJIT:
      raise RuntimeError("Not yet supported") 
    GenSweptParams(varDict)# var1Name,var1Vals)
  else: 
    if useJIT:
      runParamsFast(varDict=varDict,\
              odeName = odeName,
              name=name,dtn=deltaT,dt=dt, stim_period = stim, downsampleRate = downsampleRate)
    else:
      raise RuntimeError("antiquated") 
      runParams(runner=runner,varDict=varDict,\
              name=name,deltaT=deltaT,dt=dt, stim_period = stim, downsampleRate = downsampleRate)
  
