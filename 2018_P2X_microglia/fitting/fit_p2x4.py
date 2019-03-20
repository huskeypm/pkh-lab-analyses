"""
Demo of fitting algorithm for p2x4 fitting.
Note: the exptl data used for comparison was shifted and rescaled in an ipython notebook, so it's better to do this visually. 
"""
import matplotlib.pylab as plt 
# routlines for analyzing odes
import analyzeGotran as ao
import fittingAlgorithm as ga
fileName = "p2x4_fit.png"
testState = "I_f"

import numpy as np
scale = 1/30. * 1e-9
expdata = np.loadtxt("../validation/Toulme2012_Fig6c.csv",skiprows=4, delimiter=",")
expdata[:,1]*=scale
plt.plot(expdata[:,0],expdata[:,1])


expInterpTs = np.linspace(0.3,30,4)  # point at t=0 was pretty bad
expInterpVals = np.interp(expInterpTs,expdata[:,0],expdata[:,1])


# here I apply the offset to align with the model predictions
# also the exptl data is in seconds 
expInterpTs*= 1e3
offset = 1.55e3 # ms
expInterpTs+= offset

results = ga.run(
    odeModel = "microgliav48.ode",
    yamlVarFile = "P2X4FitVar.yaml" ,
    myVariedParam = "H2_ptxf",
    variedParamTruthVal = 0.00026  , # was 2.22e-13 in code, but bumping it up a bit
    fileName = fileName,
    timeStart = offset,
    jobDuration = 35e3,
    numRandomDraws = 15,
    numIters = 10, 
    sigmaScaleRate = 0.45,
    outputParamName = "I",
    outputParamSearcher = testState,
    outputParamMethod = "val_vs_time",
    outputParamTruthTimes=expInterpTs,
    outputParamTruthVal=expInterpVals,
    debug = True
)


