#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
__author__ = "Johan Hake (hake.dev@gmail.com)/pkh"
__date__ = "2013-03-13 -- 2014-06-12"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU LGPL Version 3.0 or later"

from modelparameters.codegeneration import latex
try:
    from scipy.optimize import root
except:
    root = None
from itertools import cycle
#import matplotlib.pyplot as plt
import numpy as np
import instant
from gotran.model.loadmodel import load_ode
from gotran.model.utils import DERIVATIVE_EXPRESSION, special_expression
from gotran.codegeneration.compilemodule import compile_module
from gotran.common.options import parameters
from gotran.common import error, warning

from modelparameters.parameterdict import *

def main(filename, params,tsteps=None):

    # Copy of default parameters
    generation = parameters.generation.copy()

    # Set body parameters
    generation.code.body.representation = params.code.body_repr
    generation.code.body.use_cse = params.code.use_cse
    generation.code.body.optimize_exprs = params.code.optimize_exprs

    # Set what code we are going to generate and not
    for what_not in ["componentwise_rhs_evaluation",
                     "forward_backward_subst",
                     "linearized_rhs_evaluation",
                     "lu_factorization",
                     "jacobian"]:
        generation.functions[what_not].generate = False

    # Always generate code for monitored expressions
    generation.functions.monitored.generate = True

    # If scipy is used to solve we generate RHS and potentially a jacobian
    if params.solver == "scipy":
        generation.functions.rhs.generate = True
        generation.functions.jacobian.generate = params.code.generate_jacobian
    else:
        generation.solvers[params.solver].generate = True

      
    # Compile executeable code from gotran ode
    ode = load_ode(filename)

    # Check for DAE
    if ode.is_dae:
        error("Can only integrate pure ODEs. {0} includes algebraic states "\
              "and is hence a DAE.".format(ode.name))
    
    # Get monitored and plot states
    plot_states = params.plot_y

    # Get x_values
    x_name = params.plot_x
    
    state_names = [state.name for state in ode.full_states]
    monitored_plot = [plot_states.pop(plot_states.index(name)) \
                      for name in plot_states[:] if name not in state_names]

    monitored = []
    for expr in sorted(ode.intermediates + ode.state_expressions):
        if expr.name not in monitored:
            monitored.append(expr.name)

    # Check valid monitored plot
    for mp in monitored_plot:
        if mp not in monitored:
            error("{} is not a state or intermediate in this ODE".format(mp))
        
    # Check x_name
    if x_name not in ["time"]+monitored+state_names:
        error("Expected plot_x to be either 'time' or one of the plotable "\
              "variables, got {}".format(x_name))

    # Logic if x_name is not 'time' as we then need to add the name to
    # either plot_states or monitored_plot
    if x_name != "time":
        if x_name in state_names:
            plot_states.append(x_name)
        else:
            monitored_plot.append(x_name)

    module = compile_module(ode, params.code.language, monitored, generation)

    parameter_values = params.parameters
    init_conditions = params.init_conditions

    if len(parameter_values) == 1 and parameter_values[0] == "":
        parameter_values = []

    if len(init_conditions) == 1 and init_conditions[0] == "":
        init_conditions = []

    if len(parameter_values) % 2 != 0:
        error("Expected an even number of values for 'parameters'")

    if len(init_conditions) % 2 != 0:
        error("Expected an even number of values for 'initial_conditions'")

    user_params = dict()
    for param_name, param_value in [(parameter_values[i*2], parameter_values[i*2+1]) \
                                    for i in range(int(len(parameter_values)/2))]:
        
        user_params[param_name] = float(param_value)

    user_ic = dict()
    for state_name, state_value in [(init_conditions[i*2], init_conditions[i*2+1]) \
                                    for i in range(int(len(init_conditions)/2))]:
        
        user_ic[state_name] = float(state_value)

    # Use scipy to integrate model
    t0 = 0.
    t1 = params.tstop
    dt = params.dt

    rhs = module.rhs
    jac = module.compute_jacobian if params.code.generate_jacobian else None
    y0 = module.init_state_values(**user_ic)
    model_params = module.init_parameter_values(**user_params)

    # Check for steady state solve
    if params.steady_state.solve and root:
        result = root(rhs, y0, args=(0., model_params,), jac=jac, \
                      method=params.steady_state.method, \
                      tol=params.steady_state.tol)

        if result.success:
            y0 = result.x
            print ("Found stead state:", ", ".join("{0}: {1:e}".format(\
                state.name, value) for value, state in zip(y0, ode.full_states)))
        else:
            warning(result.message)

    # use uniform steps 
    if type( tsteps ) is not np.ndarray:
      tsteps = np.linspace(t0, t1, t1/dt+1)

    # What solver should we use
    if params.solver == "scipy":
        try:
            from scipy.integrate import odeint
        except Exception(e):
            error("Problem importing scipy.integrate.odeint. {}".format(e))
        # OLD 
        mxstep = 10000 
        results = odeint(rhs, y0, tsteps, Dfun=jac, args=(model_params,),mxstep=mxstep)#, hmax=.03,rtol=1e-12, atol=1e-12)

    else:

        # Get generated forward method
        forward = getattr(module, "forward_"+params.solver)
        #forward = getattr(module, params.solver)
        
        results = [y0]
        states = y0.copy()
        
        # Integrate solution using generated forward method
        for ind, t in enumerate(tsteps[:-1]):

            # Step solver
            forward(states, t, dt, model_params)
            results.append(states.copy())

    # Plot results
    #if not (plot_states or monitored):
    #    return
    
    return results,module,tsteps,model_params,ode 

#if __name__ == "__main__":
def init():
    import sys, os

    body_params=parameters.generation.code.body.copy()
    
    code_params = ParameterDict(
        language = OptionParam("C", ["Python", "C"]),
        body_repr = dict.__getitem__(body_params, "representation"),
        use_cse = dict.__getitem__(body_params, "use_cse"),
        optimize_exprs = dict.__getitem__(body_params, "optimize_exprs"),
        generate_jacobian = Param(False, description="Generate and use analytic "\
                                  "jacobian when integrating."),\
        )
    
    steady_state = ParameterDict(
        solve = Param(False, description="If true scipy.optimize.root is used "\
                      "to find a steady state for a given parameters."),
        method = OptionParam("hybr", ["hybr", "lm", "broyden1", "broyden2", \
                                      "anderson", "linearmixing", "diagbroyden", \
                                      "excitingmixing", "krylov"]),
        tol = ScalarParam(1e-15, description="Tolerance for root finding algorithm."),
        )

    solver = OptionParam("scipy", ["scipy"]+\
                         list(parameters.generation.solvers.keys()),
                         description="The solver that will be used to "\
                         "integrate the ODE.")

    scipy_params = ParameterDict(mxstep=10000, hmax=.03, rtol=1e-12, atol=1e-12)
    
    params = ParameterDict(\
        scipy_params = scipy_params,
        solver = solver,
        steady_state = steady_state,
        parameters = Param([""], description="Set parameter of model"),
        init_conditions = Param([""], description="Set initial condition of model"),
        tstop = ScalarParam(100., gt=0, description="Time for stopping simulation"),\
        dt = ScalarParam(0.1, gt=0, description="Timestep for plotting."),\
        plot_y = Param(["V"], description="States or monitored to plot on the y axis."),\
        plot_x = Param("time", description="Values used for the x axis. Can be time "\
                       "and any valid plot_y variable."), 
        code = code_params)

    
    # PKH 
    #params.parse_args(usage="usage: %prog FILE [options]")
   # 
    #if len(sys.argv) < 2:
    #    error("Expected a single gotran file argument.")

    #if not os.path.isfile(sys.argv[1]):
    #    error("Expected the argument to be a file.", exception=IOError)
	 
    #file_name = "shannon_2004.ode"
    #results=main(file_name, params)

    return params 

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
def plotit(resultsName="results.txt"):
  import matplotlib.pylab as plt
  pltName = "test.png"
  print ("Printing", pltName)
  results = np.loadtxt(resultsName)
  plt.plot(results)
  plt.gcf().savefig(pltName)     

def doit(resultsName="results.txt",
         tstop=5000,ks=16.,KSRleak=5.3e-6,stim_period=1000, # ms
         file_name = "shannon_2004.ode",
         stateName = "Cai"
         ):
  print ("Using ode model ", file_name )
  print ("Writing results to ", resultsName)
  params = init()
  params.tstop = tstop 

#  file_name = "shannon_2004.ode"
  params.parameters = ['ks',ks,'KSRleak',KSRleak,'stim_period',stim_period]
  results,module,tsteps,model_params,ode=main(file_name, params)

  # save cai 
  idx =  module.state_indices( stateName ) 
  np.savetxt(resultsName,results[:,idx])


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -doit/-plotit" % (scriptName)
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

  # Loops over each argument iln the command line 
  file_name = "shannon_2004.ode"
  resultsName = "results.txt" 
  stateName = "Cai"
  stim_period = 1000. # ms 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-doit"):
      doit(file_name = file_name,resultsName=resultsName,
           stateName=stateName, stim_period = stim_period)
      quit()
    if(arg=="-2"):
      doit(resultsName="mod.txt",ks=14,KSRleak=5.2e-6)
      quit()
    if(arg=="-plotit"):
      plotit()
      quit() 
    if(arg=="-stim_period"):
      stim_period = np.float(sys.argv[i+1])
    if(arg=="-file_name"):
      file_name = sys.argv[i+1]
    if(arg=="-resultsName"):
      resultsName= sys.argv[i+1]
    if(arg=="-stateName"):
      stateName= sys.argv[i+1]
  





  raise RuntimeError("Arguments not understood")




