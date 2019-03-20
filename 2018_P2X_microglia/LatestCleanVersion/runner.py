from separate import * # This imports separate fluxes
import shannon_2004_hack as model
# NOTE: this hack has the cut-and-pasted parameters, as described in 
# 150220_gotranbug 

## MONITORS
# New gotran 
#monitors = ["fCa_SL","fCa_jct","i_NaCa","j_rel_SR","j_pump_SR","i_Stim"]
#units = ["unk","unk","unk","unk","unk","unk"]
#idxMonitors = np.zeros(np.shape(monitors)[0],dtype=int)
#for i,monitor in enumerate(monitors):
#  idxMonitors[i] = model.monitor_indices(monitor)


## STATE VAR
Cai_idx =  model.state_indices( "Cai" )
Ca_SR_idx =  model.state_indices( "Ca_SR" )
V_idx =  model.state_indices( "V" )


## PARAMS
# stim_period=(121
#stim_period_pIdx=121
#V_max_Jpump_pIdx = 71 # SERCA
#V_max_pIdx = 45 # NCX

# old gotran 
#stim_period_pIdx =  model.parameter_indices( "stim_period" )
#V_max_pIdx =  model.parameter_indices( "V_max" )

# new gotran
stim_period_pIdx =  model.parameter_indices( "stim_period" )
#V_max_pIdx =  model.parameter_indices( "V_max" )



## Misc
mM_to_uM = 1e3

## Monitors 
#huskeypm@huskeypm-ubuntu12:~/sources/wholecell$ grep monitor shannon_2004.ode 
#monitor(fCa_SL) 
#monitor(fCa_jct) 
#monitor(i_NaCa)   
#monitor(j_rel_SR)
#monitor(j_pump_SR)
#monitor(i_Stim)
#fCa_SL  =0
#fCa_jct =1
#i_NaCa  =2 
#j_rel_SR=3
#j_pump_SR=4
#i_CaL=5
#i_Stim = 6
#totMonitors = 7

from scipy.integrate import odeint
import numpy as np

# get monitored values 
def monitorstepper(model,states,pi,tsteps):
  dtt= tsteps[1]-tsteps[0]
  tstepsm1 = tsteps[1::] 

  # get len
  si = states[0,:]
  jis = model.monitor(states[0,:], 0, pi) 
  totMonitors = np.shape(jis)[0]
  jall = np.zeros((np.shape(tstepsm1)[0],totMonitors)) 
  idxNCX = model.monitor_indices("i_NaCa")

  #print np.shape(jis)
  jSums = np.zeros(totMonitors)
  for i,t in enumerate(tstepsm1):
    # get current state
    si = states[i,:]
    # extract ALL monitored fluxes 
    jis = model.monitor(si, t, pi) 
    #print np.shape(jis)
    jall[i,] = jis
    
    # sum fluxes 
    jSums += jis*dtt   
    
    
  return (tstepsm1,jall)  

def init():
  # old gotran 
  #s=model.init_values()
  #p=model.default_parameters()

  # new gotran 
  s=model.init_state_values()
  p=model.init_parameter_values()
  t=0; #dt=1000; dtn=5;

  model.s = s; model.t =t; model.p = p


# run simulation 
# dtn <1, otherwise action potential isn't correct
def runner(dt=1000,dtn=1.,\
           stim_period=1000,\
  #         V_max_Jpump = 0.0053114, # SERCA, [mM/ms]
  #         V_max = 9, # NCX, [uA/uF]
           mxsteps=10000):

  s = model.s; t =model.t; p = model.p
  
  # Important variables, states
  p[stim_period_pIdx]=stim_period
  #p[V_max_Jpump_pIdx]= V_max_Jpump 
  #p[V_max_pIdx] = V_max


  
  # Basic run and grab outcomes
  tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
  
  #states = odeint(model.rhs,s,tsteps,(p,),mxstep=mxsteps)
  #print "USING SUPPER!" 
  # mxsteps 
  states = odeint(model.rhs,s,tsteps,(p,),mxstep=mxsteps,hmax=.03,rtol=1e-12, atol=1e-12)
  #V_idx =  model.state_indices( "V" )
  #plt.plot(tsteps,mM_to_uM*states[:,V_idx],label=Cai_idx)
  
  # get monitored variables
  (ts,js)=monitorstepper(model,states,np.copy(p),tsteps)
  ts=tsteps[1:]; #js=0
  
  return (p,states,ts,js)


def plotting(p,states,ts,js,case="default"):
  import matplotlib.pyplot as plt
  tsteps = np.zeros(np.shape(ts)[0]+1)
  tsteps[1:]=ts
  state_indices = model.state_indices    


  ## Ca transients 
  plt.figure()
  plt.subplot(2,1,1)
  plt.plot(tsteps,mM_to_uM*states[:,Cai_idx],label=Cai_idx)
  plt.title("Cytosolic Ca")
  plt.ylabel("Ca [uM]")
  plt.xlabel("t [ms]") 
  
  plt.subplot(2,1,2)
  plt.plot(tsteps,mM_to_uM*states[:,Ca_SR_idx],label=Ca_SR_idx)
  plt.title("SR Ca")
  plt.ylabel("Ca [uM]")
  plt.xlabel("t [ms]") 
  
  plt.subplots_adjust(hspace=0.5) 
  plt.gcf().savefig(case+"_calcium.png",dpi=300)
  
  ## voltage 
  plt.figure()
  plt.plot(tsteps,states[:,V_idx],label=V_idx)
  plt.ylabel("V [mV]")
  plt.xlabel("t [ms]") 
  plt.gcf().savefig(case+"_potential.png",dpi=300)

  return 
  
  ## fluxes
  # BROKEN (ts,js)=monitorstepper(model,states,np.copy(p),tsteps)

  #
  plt.figure(figsize=(10,5))
  plt.subplot(1,2,1)
  plt.plot(ts,js[:,i_CaL_idx],label="i_CaL")               
  plt.legend(loc=0)
  
  plt.subplot(1,2,2)
  plt.plot(ts,js[:,i_NaCa_idx],label="i_NaCa_idx")
  plt.legend(loc=0)

  #
  plt.figure(figsize=(10,5))
  plt.subplot(1,2,1)
  plt.plot(ts,js[:,j_rel_SR_idx],label="j_rel_SR_idx")
  plt.legend(loc=0)
  
  plt.subplot(1,2,2)
  plt.plot(ts,js[:,j_pump_SR_idx],label="j_pump_SR_idx")
  plt.legend(loc=0)

  # w
  plt.figure(figsize=(10,5))
  if 1:
    plt.subplot(1,2,1)
    plt.plot(ts,js[:,fCa_SL_idx],label="fCa_SL_idx")
    plt.plot(ts,js[:,fCa_jct_idx],label="fCa_jct_idx")
    plt.legend(loc=0)
  
    plt.subplot(1,2,2)
    plt.plot(ts,js[:,i_Stim_idx],label="i_Stim_idx")
    plt.legend(loc=0)

  #plt.gcf().savefig(case+"_fluxes.png",dpi=300)

  # returning only 
  return js[:,j_rel_SR_idx]

# <codecell>



#!/usr/bin/env python
import sys
#
# Revisions
#       10.08.10 inception
#

def doit(caseName="default"):
  init()
  (p,si,tsi,jsi) = runner(dt=2000,dtn=1.0,stim_period=901,mxsteps=10000)
  plotting(p,si,tsi,jsi,case=caseName)
  raise RuntimeError("Not actually a validation yet")

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

if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  caseName = "default"
  for i,arg in enumerate(sys.argv):
    if(arg=='caseName'):
      caseName=sys.argv[i+1]
    if(arg=="-validation"):
      #arg1=sys.argv[i+1] 
      doit(caseName=caseName)
      quit()
  





  raise RuntimeError("Arguments not understood")




