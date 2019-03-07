"""
Validation script for fitting algorithm. Mostly checks for consistency 
"""
import matplotlib.pylab as plt 
# routlines for analyzing odes
import analyzeGotran as aG
import numpy as np 


import fittingAlgorithm as ga


# ### Fitting a time-dependent profile (vector) 
# A more complicated example in which a time-dependent profile is compared.
# Below we first load the reference data 
# 
# <font color=red>
# WARNING: I'm rescaling the Toulme data by a factor of 40, since their transgenic model behaves as if it is activated and we'd like values more apporpriate for resting glia. The 1e-9 factor is for converting into uC/ms <--- PLease verify
# </font>
# 

# In[6]:

scale = 1/30. * 1e-9
expdata = np.loadtxt("../validation/Toulme2012_Fig6c.csv",skiprows=4, delimiter=",")
expdata[:,1]*=scale


# We chose an offset to account for any difference between the exp. data and when the channel fires in our data. This took some playing around. Easiest to debug this when looking at the best fit graph at the end of the calculation. Also, I interpolated the experimental over about 10 s. Need to also play around with interpolation points too so they fairly approximate the complete data 

# In[7]:


expInterpTs = np.linspace(0.3,10,4)  # point at t=0 was pretty bad
expInterpVals = np.interp(expInterpTs,expdata[:,0],expdata[:,1])
plt.plot(expdata[:,0],expdata[:,1])
plt.scatter(expInterpTs,expInterpVals)


# here I apply the offset to align with the model predictions
# also the exptl data is in seconds 
expInterpTs*= 1e3
offset = 1.55e3 # ms
expInterpTs+= offset

#print expInterpTs

def singleVarTest():
  testState = "I_ptxf"
  results = ga.run(
      odeModel = "microgliav48.ode",
      yamlVarFile = "P2X4FitVar.yaml" ,
      myVariedParam = "G12_ptxf",
      variedParamTruthVal = 4.22e-13  , # was 2.22e-13 in code, but bumping it up a bit
      timeStart = offset,
      jobDuration = 20e3,
      numRandomDraws = 10,
      numIters = 10,
      sigmaScaleRate = 0.15,
      outputParamName = "I",
      outputParamSearcher = testState,
      outputParamMethod = "val_vs_time",
      outputParamTruthTimes=expInterpTs,
      outputParamTruthVal=expInterpVals,
      debug = True
  )
  
  
  
  ##
  ## Assert to check that responses still ok 
  ##
  
  data = results['data']
  stateLabel = "I_ptxf"
  subData = aG.GetData(data,stateLabel)
  ts = subData.t
  vals = subData.valsIdx
  
  
  predInterpVals= np.interp(expInterpTs*1e-3,ts,vals)
  
  err = np.linalg.norm(predInterpVals -expInterpVals)
  prevErr = 6.78742768088e-10 * 1.1
  print err,prevErr
  assert(err < prevErr), "Something changed!"
  print "PASS!"

def multiVarTest(iters=10,numRandomDraws=24): 
  testState = "I_f"
  # Best values from before 
  #    k1_ptxf = ScalarParam(1.0e-4, unit="ms**-1"),
  #    k2_ptxf = ScalarParam(2610, unit="(M*ms)**-1"),
  #    k3_ptxf = ScalarParam(3.915, unit="ms**-1"),
  #    k4_ptxf = ScalarParam(11537.9, unit="(M*ms)**-1"),
  #    k5_ptxf = ScalarParam(0.000652, unit="ms**-1"), #3.0e1
  #    k6_ptxf = ScalarParam(3461.38, unit="(M*ms)**-1"),
  #    H1_ptxf = ScalarParam(1.67e-5, unit="ms**-1"),
  #    H2_ptxf = ScalarParam(0.00026, unit="ms**-1"),
  #    H3_ptxf = ScalarParam(1.56e-3, unit="ms**-1"),
  #    H4_ptxf = ScalarParam(0.0237e-3, unit="ms**-1"),
  #    H5_ptxf = ScalarParam(0.0013, unit="ms**-1"), #8e-3
  #    H6_ptxf = ScalarParam(0.00018, unit="ms**-1")
  stddev = 0.2
  variedParamDict = {
      # paramDict[myVariedParam] = [variedParamTruthVal, 0.2] # for log normal
      'k1_ptxf': [1.0e-4, stddev], 
      'k2_ptxf': [2610, stddev],
      'k3_ptxf': [3.915, stddev], 
      'k4_ptxf': [11537.9, stddev],
      'k5_ptxf': [0.000652, stddev], 
      'k6_ptxf': [3461.38, stddev],
      'H1_ptxf': [1.67e-5, stddev], 
      'H2_ptxf': [0.00026, stddev], 
      'H3_ptxf': [1.56e-3, stddev], 
      'H4_ptxf': [0.0237e-3, stddev], 
      'H5_ptxf': [0.0013, stddev], 
      'H6_ptxf': [0.00018, stddev]
  }
  results = ga.run(
      odeModel = "microgliav48.ode",
      yamlVarFile = "P2X4FitVar.yaml" ,
      variedParamDict = variedParamDict,
      timeStart = offset,
      jobDuration = 35e3,
      numRandomDraws = numRandomDraws,
      numIters = iters,
      sigmaScaleRate = 0.05,
      outputParamName = "I",
      outputParamSearcher = testState,
      outputParamMethod = "val_vs_time",
      outputParamTruthTimes=expInterpTs,
      outputParamTruthVal=expInterpVals,
      debug = True
  )

  previousFitness = 0.161 + 0.01
  assert(results['bestFitness'] < previousFitness), "SOmethign changed" 

# tests that non-uniform time steps give correct answers
import runModel as rs
def nonunifTime():
  ## Full 
  odeName="microgliav55.ode"
  yamlVarFile="nfatFitVar.yaml"
  dt=0.1
  #dtn =900e3  # job duration 
  dtn = 20e3  # job duration 
  stateLabel = "NFATNn"
  #stateLabel = "I_ptxf"
  fixedParamDict = aG.YamlToParamDict(yamlVarFile)
  
  outName = "full"
  returnDict = dict() # return results vector
  if 1: 
    rs.runParamsFast(odeName=odeName,
                       name=outName,
                       varDict = fixedParamDict,
                       dt=0.1,dtn=dtn,
                       generateJacobian=True,
                       returnDict=returnDict)
    data = aG.readPickle(outName+".pkl")
  
  ## non-unif
  tsteps = rs.TimeSteps(dtFine=dt,dtnFine=3e3,dtCoarse=1.0e3,dtnCoarse=dtn)
  outName = "nonunif"
  returnDictNU = dict() # return results vector
  rs.runParamsFast(odeName=odeName,
                       name=outName,
                       varDict = fixedParamDict,
                       tsteps = tsteps,
                       generateJacobian = True,
                       returnDict=returnDictNU)
  dataNU = aG.readPickle(outName+".pkl")
  
  ## Compare 
  
  subDataFull = aG.GetData(data,stateLabel)
  plt.figure()
  plt.plot(subDataFull.t,subDataFull.valsIdx, label="full, %d points"%np.shape(subDataFull.t)[0])
  #print 
  
  subData = aG.GetData(dataNU,stateLabel)
  valsFull = np.interp(subData.t,subDataFull.t,subDataFull.valsIdx)
  err = np.linalg.norm(subData.valsIdx - valsFull)

  plt.plot(subData.t,subData.valsIdx,'k.', label="nonunif, %d points\nerr=%e"%(np.shape(subData.t)[0],err))
  plt.legend(loc=0)
  plt.title(stateLabel)
  plt.xlabel("time [s]") 
  plt.gcf().savefig("comp.png") 
  
  assert(err < 1e-10), "Non uniform time step broke!"



def test():
  nonunifTime()



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
def validate():
  # Run tests
  singleVarTest()
  multiVarTest()
  nonunifTime()


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
      validate()
      quit()
    if(arg=="-multi"):
      multiVarTest(iters=30,numRandomDraws=64)
      quit()
    if(arg=="-test"):   
      test()
      quit()
  





  raise RuntimeError("Arguments not understood")




