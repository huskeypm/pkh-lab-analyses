# Function for integreating ODE and displaying results
import scipy.integrate
from scipy.integrate import odeint
import numpy as np
import matplotlib.pylab as plt
import math
from math import exp
from scipy.interpolate import spline
import pickle as pk


def compareProfiles(lit_data,model_data): # indpnt = independent variable, dpnt = dependent variable
    lit_xset, lit_yset = lit_data    # literature data
    model_xset, model_yset = model_data  # modelled data

    est_yset = np.interp(lit_xset,model_xset,model_yset)

    N = np.shape(est_yset)[0]
    iters = np.arange(N)

    #err2 = 0

    #for i in iters:
    #    err2 = err2 + (lit_yset[i] - est_yset[i])**2
    err2 = np.sum( (lit_yset - est_yset)**2) 

    err2 = err2*100
    return err2
