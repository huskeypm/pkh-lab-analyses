"""
Main program for running MSA approximation for chemical potentials of hard charged spheres in a binding site 
"""

import numpy as np 
## Constants and conversions 
N = 6.022e23 # #/L
L_to_mL = 1e-3 
mL_to_nm3 = (1e-7)**3 # 1/cm^3 * (cm/nm)**3
M_to_N_o_nm3  = L_to_mL * mL_to_nm3 * N 
epsilonFilter = 63.5 # from Nonner
kT_mV = 25.6 # mV * e0
e0 = 1/kT_mV # 1/mV
nOxy = 8.

# In[2]:

import MSAtools as msa

def doit(mode="nonner",
         iters = 1000,
         nOxy=8):
  ## our system's config
  #mode="ours"
  if mode=="ours":
    conc_M = np.array([1e-20, (150.e-3+ 0.2e-6 + 4e-3), 150e-3, 0.1e-6, 2e-3])  # [M] order is  O Cl, K, Mg, Ca (bath) 
    zs = np.array([-0.5, -1.0, 1.0, 2.0, 2.0])
    Ns = np.array([6., 1., 1., 1., 1.]) # Filter
    V_i = 0.16464 #3.1 angstrom radius #0.1438
    sigmas = np.array([0.278, 0.362, 0.276, 0.200, 0.250])
  
  elif mode=="nonner":
    msa.idxOxy = 0
    #msa.idxNa = 2
    conc_M = np.array([1e-20, (100.0e-3), 100.0e-3, 1.0e-10])  # [M] order is  O Cl, Na, Mg, Ca (bath) 
    zs = np.array([-0.5, -1.0, 1.0, 2.0])
    Ns = np.array([8, 0, 2, 1]) # Filter
    V_i = 0.375 #3.1 angstrom radius #0.1438
    sigmas = np.array([0.278, 0.362, 0.204, 0.200])

  ## run solver 
  rhoFilter,donnanPotential = SolveMSAEquations(conc_M,zs,Ns,V_i,sigmas)

def SolveMSAEquations(conc_M,zs,Ns,V_i,sigmas,
  muiexsPrev = 0.,
  psiPrev = 0.,
  alpha = 1e-3,  # convergence rate (faster values accelerate convergence) 
  maxiters = 1e5, # max iteration before loop leaves in dispair 
  verbose=False): 
  
  ## Convert concs/numbers into densities 
  conc_M_to_N_o_nm3 = conc_M * M_to_N_o_nm3 # [#/nm3]
  print conc_M_to_N_o_nm3
  rhoisBath = conc_M_to_N_o_nm3
  rhoisFilter = np.array(Ns / V_i)
  
  
  
  # Compute debye length 
  #kappainv = 0.304 / np.sqrt(conc_M[idxNa]) # Israelach
  #print "Analy kinv [nm] %f (nacl only)" %(kappainv) # for 100 mM
  kappa = msa.CalcKappa(conc_M) # np.sqrt(4*np.pi *msa.lambda_b * (2*6.022e26* conc_M[idxNa]))
   
  print "Est kinv [nm] %f " % (1e9/kappa)
  #print kappa
  
  ## compute gamma for bath 
  GammaBath = kappa/2*(1/msa.m_to_nm)   
  muexBath_kT = msa.CalcMuexs(GammaBath,zs)/msa.kT_kcalmol
  # Validated (with my notes)
  #print "mu_Na,bath [kT] " , muexBath_kT[idxNa]
  print "muexBath_kT ", muexBath_kT
  
  
  ## Calculate Bjerrum, filter chem potential 
  
  #GammaFilter = 2.11 # [nm] 
  # verified
  #lambdaBFilter = msa.CalcBjerrum(epsilon=epsilonFilter)
  #print lambdaBFilter*msa.m_to_nm 
  #muFilter =  msa.CalcMuexs(GammaFilter,zs,\
  #                          lambda_b = lambdaBFilter)/msa.kT_kcalmol
  
  #print "muFilter_Na [kT] ", muFilter[msa.idxNa]
  
  #psi = -170 # [mV] from Nonner for 0.1 M NaCl
  #elecPotentContrib = zs[idxNa]*e0*psi
  #
  ## consistent with notes
  #print "z e0 psi [kT] ", elecPotentContrib
  
  
  ## Calculate Hardsphere components 
  def xiStuff(rhos,sigmas): 
    xi_3_= msa.xis(rhos, sigmas, 3)
    delta_= 1. - xi_3_
    xi_1_= msa.xis(rhos, sigmas, 1)
    xi_2_= msa.xis(rhos, sigmas, 2)
    xi_0_ = msa.xis(rhos, sigmas, 0)
    #print "xi_n (n= 0 to 3) =", xi_0_, xi_1_, xi_2_, xi_3_; 
    #print "delta =", delta_
    return xi_0_, xi_1_, xi_2_, xi_3_, delta_

  xi_0_filter, xi_1_filter, xi_2_filter, xi_3_filter, delta_filter = xiStuff(rhoisFilter,sigmas)
  print "xi_n (n= 0 to 3) =", xi_0_filter, xi_1_filter, xi_2_filter, xi_3_filter
  print "delta =", delta_filter
  xi_0_bath, xi_1_bath, xi_2_bath, xi_3_bath, delta_bath = xiStuff(rhoisBath,sigmas)
  print "xi_n (n= 0 to 3) =", xi_0_bath, xi_1_bath, xi_2_bath, xi_3_bath
  print "delta =", delta_bath
  
  
  ## Compute deltaphi components 
  HSphifilter = msa.deltaphiHS(xi_0_filter, xi_1_filter, xi_2_filter,xi_3_filter,delta_filter)
  #print HSphifilter
  
  HSphibath = msa.deltaphiHS(xi_0_bath, xi_1_bath, xi_2_bath,xi_3_bath,delta_bath)
  
  
  ## Compute HS contribution to chem potential for filter and bath 
  mu_HS_filter = msa.HSmuis(sigmas, delta_filter, xi_0_filter, xi_1_filter, xi_2_filter,xi_3_filter, HSphifilter)
  mu_HS_bath = msa.HSmuis(sigmas, delta_bath, xi_0_bath, xi_1_bath, xi_2_bath, xi_3_bath, HSphibath)
  
  #print mu_HS_filter
  #print muexBath_kT
  muexBath_kT += mu_HS_bath
  print "muexBath_kT ", muexBath_kT
  
  
  nIons = np.shape(conc_M)
  Vs = np.ones(nIons)*V_i
  lambdaBFilter = msa.CalcBjerrum(epsilon=epsilonFilter)
  import scipy.optimize
  
  ## Main LOOP
  print "MAIN"
  iters = 0 
  #for j in np.arange(maxiters):    
  tol = 1.0e-8
  psi = -10.
  psiDiff  = psiPrev - psi

  while(abs(psiDiff) > tol):

    ##  Get Updated Rhois              
    
    rhoisFilter = msa.CalcRhoFilter(muexBath_kT,muiexsPrev,psiPrev,zs,rhoisBath)
    rhoisFilter[msa.idxOxy] = nOxy/V_i
    if verbose: 
      print "p [M] (Before  oxy correct): ", rhoisFilter
      print "p [M] (After oxy correct): ", rhoisFilter
  #if 0:            
  #            rhoisFilter = CalcRhoiFilter(muBath,muiexsPrev,psiPrev,zis,rhoisBath)            
  #            rhoOxy = nOxy / V_body
  #            rhoisFilter[iOxy] = rhoOxy
  #            print "pFilter",rhoisFilter
  #            print "WARNING: overriding updated rho_Oxy w fixed value"
          
      # Vectorize
      # vfunc = np.vectorize(CalcRhoiFilter)
      # rhoisFilter = vfunc(muisBath,...)
  
      ## Updated Muis via MSA 
          
    #x=[GammaPrev] # from filter 
  
      # Get updated gamma
      # Passing rho as N/V, therefor using Ns = rho, Vs = 1
    #GammaFilter = MSAeqn(x,rhoisFilter,np.ones(nIons))
    #kappa = sqrt(4 * np.pi * lambda_b * M_to_nm * sum(rhoisFilter*zis*zis))
    #Gamma = kappa/2. 
    GammaFilter = msa.CalcKappa(rhoisFilter)/2.  # use Gamma based on filter concs
    GammaFilter/= msa.m_to_nm
    #print GammaFilter
    xopt = scipy.optimize.fmin(func=msa.MSAeqn,x0=GammaFilter,disp=False,\
                               args=(rhoisFilter,Vs,zs,sigmas,lambdaBFilter))
    GammaFilter = xopt[0]
  
    #print "GammaF %f" % GammaFilter
      # get updated muis    
    #muiexs = np.zeros(nIons)
    #print zis
    
    
    mu_ES = msa.CalcMuexs(GammaFilter,zs,lambdaBFilter)
    mu_HS = msa.CalcHS(rhoisFilter,sigmas)    
    muiexs = mu_ES + mu_HS
    if verbose: 
      print "muiexs " , muiexs      # relax muis
      #### s.b. eqn 31
    muiexs = (alpha*muiexs + muiexsPrev)/(1.+ alpha) 
    if verbose: 
      print "rescaled muiexs " , muiexs
    
      ## Update psi
    psi = msa.UpdatePsi(rhoisFilter,zs,psiPrev)
  
    psiDiff = np.abs(psiPrev-psi)/psi
    if verbose: 
      print "psi %f/psiDiff %f"% (psi,psiDiff)
      print "rhoisFilter ", rhoisFilter
      print "======="
      ## Update prev
    psiPrev = psi
    muiexsPrev = muiexs
    iters += 1
    if iters > maxiters:
      print "Your function broke!"
      break 
  # muiexs - chem potential of filter 
  # psi donnan potential 
    ## END LOOP 
  return muiexs,psi,mu_ES,mu_HS,rhoisFilter
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
      #arg1=sys.argv[i+1] 
      doit(mode="nonner") # arg1)
      quit()
    if(arg=="-run"):
      doit(mode="ours") 
      quit()
  





  raise RuntimeError("Arguments not understood")




