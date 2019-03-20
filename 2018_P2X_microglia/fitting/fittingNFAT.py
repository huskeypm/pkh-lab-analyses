import matplotlib.pylab as plt 
# routlines for analyzing odes
import analyzeGotran as aG
import numpy as np
import fittingAlgorithm as ga
import numpy as np 
import runModel as rm
# #### Fitting kinetic data
# Will assume that for an equilibrium reaction that is far from equilibrium
# A + B <--> AB
# that the contributions of kf to AB' are most significant and AB' is proportional to kf 
# THerefore, can determine kf by measuring initial slope 
# 

# Here we run everything. The best fit parameter is contained in results
# results['bestFitParam']

# In[111]:


yamlVarFile = "nfatFitVar.yaml" 
odeModel = "microgliav55.ode"
jobDuration = 900e3  # [ms]
dtnFine = 3e3        # Fine time steps from 0..dtnFine [ms]
dtFine=0.1           # Fine dt [ms]
tsteps = rm.TimeSteps(dtFine=dtFine,dtnFine=dtnFine,dtCoarse=1e3,dtnCoarse=jobDuration)


# In[112]:


# First run 'original' case

outName = "orig"
fixedParamDict = aG.YamlToParamDict(yamlVarFile)
if 1:
 rm.runParamsFast(odeName=odeModel,
                     name=outName,
                     varDict = fixedParamDict,
                     tsteps = tsteps)
dataOrig = aG.readPickle(outName+".pkl")



# now run GA


import runModel as rm
stddev = 0.02  
variedParamDict = {
    # paramDict[myVariedParam] = [variedParamTruthVal, 0.2] # for log normal
    'kf1_NFAT': [9.0e-11,stddev],  
    'kr1_NFAT': [0.53e-3, stddev],
    'kf2_NFAT': [1.44e-3,stddev],
    'kf3_NFAT': [3.62e-4,stddev], 
    'kr3_NFAT': [9.71e-5,stddev], 
    'kf4_NFAT': [4.45e-4,stddev], 
    'KmN_NFAT': [535,stddev] 
}


results = ga.run(
    odeModel = odeModel,
    yamlVarFile=yamlVarFile,
    variedParamDict = variedParamDict,
    timeStart = dtnFine+dtFine, # [ms] 
    tsteps = tsteps,
    numRandomDraws = 20,
    numIters = 20, 
    sigmaScaleRate = 0.05,
    outputParamName = "deNFAT nuc.",
    outputParamSearcher = "NFATNn",
    outputParamMethod = "tnorm_maxdydt",
    outputParamTruthVal=1.11e-6,  # 1/900e3
    debug = True
)


# my guesses for nomralized nfat expt data
expt=np.array([0,2,15])
expv=np.array([0,0.25,.95])


#plt.scatter(expInterpTs, expInterpVals)
plt.figure()
key="NFATNn"
data = results['data'] # from GA 
subData = aG.GetData(data,key)
vals = subData.valsIdx - np.min(subData.valsIdx)
#vals /= np.max(vals)
plt.plot(subData.t/(1e3*60),vals, label=key)
(1e3*60)

subData = aG.GetData(dataOrig,key)
vals = subData.valsIdx - np.min(subData.valsIdx)
origMax = np.max(vals)
#vals /= np.max(vals)
plt.plot(subData.t/(1e3*60),vals, label=key+" orig")

plt.scatter(expt,expv*origMax,label="exp")
plt.legend(loc=0)
plt.gcf().savefig("result.png") 
