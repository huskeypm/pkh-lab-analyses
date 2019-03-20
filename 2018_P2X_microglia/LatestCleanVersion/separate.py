## TEMPLATE 
#$ grep monitor shannon_2004.ode | perl -ane '$n++; chomp $_; ($nm)=($_=~/\((.*)\)/);printf("${nm}_idx=%d # $_\n",$n-1)'
# shannon PKH:
#Vol_SR_idx=0 # monitor(Vol_SR)  
#Vol_myo_idx=1 # monitor(Vol_myo) 
#fCa_SL_idx=2 # monitor(fCa_SL) 
#fCa_jct_idx=3 # monitor(fCa_jct) 
#i_NaCa_idx=4 # monitor(i_NaCa)   
#j_rel_SR_idx=5 # monitor(j_rel_SR)
#j_pump_SR_idx=6 # monitor(j_pump_SR)
#J_Ca_SL_myo_idx=7 # monitor(J_Ca_SL_myo)
#dCa_TroponinC_idx=8 # monitor(dCa_TroponinC)
#i_Stim_idx=9 # monitor(i_Stim)   
#dCa_cytosol_tot_bound_idx = 10
#totMonitored = 11
# shannon_PKH
fCa_SL_idx=0 # monitor(fCa_SL) 
fCa_jct_idx=1 # monitor(fCa_jct) 
i_NaCa_idx=2 # monitor(i_NaCa)   
j_rel_SR_idx=3 # monitor(j_rel_SR)
j_pump_SR_idx=4 # monitor(j_pump_SR)
i_CaL_idx=5 # monitor(i_CaL)   
i_Stim_idx=6 # monitor(i_Stim)   
totMonitored = 7


# This function computes the effective substrate flux 
# over a given time interval. 
# The results contain the concentration flux, the average flux due to the
# individual 'J' terms in the ode mode, as well as the averaged 'J' terms
# themselves
# not sure if this is the correct/efficient way of doing this
# runs iterator, grabs flux values at each time point 

from scipy.integrate import odeint
import numpy as np
class empty:pass 

print "WARNING: ode idx is hardcoded" 
Cai_ode_idx = 37 


## THIS IS A TEMPLATE 
## Validated that Jvol + Jbound = dca reported in shannon model 
def separateFluxes(model,s,pi,tsteps,warn=False):

  # run model to get states vs time  
  states = odeint(model.rhs,s,tsteps,(pi,))
  dt = tsteps[-1] - tsteps[0] # TODO check 

  # namely, need to determine which fluxes are going into the PDE domain 
  fluxes = empty()

  
  jSums= np.zeros(totMonitored) # total num of monitored fluxes 
  # TODO genealize 
  rhsSums= np.zeros(39) # total num of monitored fluxes 

 
  ## DONT TOUCH 
  # TODO put in its own function 
  ## grab js from each time step 
  # Idea here is that 
  # [c(Tf) - c(T0)]/(Tf-T0) = \int_T0^Tf j
  # --> delC ~ delT Sum[ j*dt] 
  # get dt between time steps
  # ignore first time step 
  dtt= tsteps[1]-tsteps[0]
  tstepsm1 = tsteps[1::]
  for i,t in enumerate(tstepsm1):
    # get current state
    si = states[i,:]
    # extract minitored fluxes 
    jis = model.monitor(si, t, pi)

    ## debugging
    rhs = model.rhs(si, t, pi)
    #print "WARING: BREAKING"; break 
    

    ## sum fluxes 
    jSums += jis*dtt
    rhsSums += rhs*dtt
    #print "m3/v",m3[7:9]/V_cyt

  ## get average jAvg over the entire time interval 
  dT =tsteps[-1] - tsteps[0]
  jAvg = jSums/dT
  rhsAvg = rhsSums/dT

  ## determine concentration flux from first/last state  
  f = states[-1,]
  s = states[0,]
  #dc= f[ATP_idx] - states[0,ATP_idx]
  dc= f - s
  fluxes.dcdt = dc/dT


  ## END LEAVE MEE 
  #print "WARNINGL: HACK WRONG!!!!!" 
  #jAvg = jis


  ## TEMPLATE
  # from ODE: 
  # dCai_dt = -j_pump_SR*Vol_SR/Vol_myo + J_Ca_SL_myo/Vol_myo - one*dCa_cytosol_tot_bound
  dCai_dt = rhs[Cai_ode_idx]
  fluxes.jBoundary = np.zeros(np.shape(s)[0]) 
  Vol_myo = jis[Vol_myo_idx] # unfortunate abuse of variable names  
  Vol_SR  = jis[Vol_SR_idx] # unfortunate abuse of variable names  
  fluxes.jBoundary[Cai_ode_idx] = jAvg[J_Ca_SL_myo_idx]/Vol_myo
  #print jAvg[J_Ca_SL_myo_idx]
  # CORRECT print "SL ", fluxes.jBoundary[Cai_ode_idx]
   

  fluxes.jVol = np.zeros(np.shape(s)[0]) 
  # this is quirkly but I think we need to substract SR component, but 
  # add in TnC component (so that we are left with just the contribution 
  # from the SL (see SB expression for dCai) 
  fluxes.jSR_novolscale = jAvg[j_pump_SR_idx]
  fluxes.jSR = -jAvg[j_pump_SR_idx]*Vol_SR/Vol_myo
  #print jAvg[j_pump_SR_idx]
  fluxes.jVol[Cai_ode_idx] = fluxes.jSR
  fluxes.jVol[Cai_ode_idx] += -jAvg[dCa_cytosol_tot_bound_idx]
  # CORRECT print "pump ", fluxes.jVol[Cai_ode_idx]

  # WARNING: not the way to do this, since i want to do buffering via TnC through PDE, not here
  tot = fluxes.jBoundary[Cai_ode_idx] + fluxes.jVol[Cai_ode_idx]
  fluxes.jVol[Cai_ode_idx] += 0 # TODO add TnC component here
  #print "here", dCai_dt, jAvg[J_Ca_SL_myo_idx], Vol_myo
  # these agree 
  #print "dCai_dt (rhs) %f vs tot %f ", (rhsAvg[Cai_ode_idx], tot)
  #print fluxes.jBoundary[Cai_ode_idx] 
  #print fluxes.jVol[Cai_ode_idx]

  
  #print Vol_myo
  #print jAvg[J_Ca_SL_myo_idx]
  #print jAvg[j_pump_SR_idx]
   
  
  return (states,fluxes)



