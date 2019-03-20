## derive Inesi model for P4
# baserd on  Kinetic and Equilibrium Characterization of an
# Energy-Transducing Enzyme and Its Partial Reactions

import numpy as np
import sympy as sp
M_to_uM = 1e6

## compare model w Fig 8 
#parameters are from pg 167
#EpCa2eqs = EpCa2eq.subs({\
#  Et: 1.,     # assume Et is fixed?? 
#  k1: 4.25e7, # 1/M sec (Af)
#  kn1: 450,   # 1/sec
#  k2: 15,     # 1/sec (Bf) 
#  kn2: 33,    # 1/sec
#  k3: 1e8,    # 1/M sec (Cf) 
#  kn3: 16     # 1/sec
#  })
#EpCa2func = sp.solve(EpCa2eqs,EpCa2)[0]

# Inesi
class paramsInesi:
    k1= 4.25e7 # 1/M sec (Af)
    kn1= 450.  # 1/sec
    k2= 15.    # 1/sec (Bf) 
    kn2= 33.   # 1/sec
    k3= 1e8    # 1/M sec (Cf) 
    kn3= 16.   # 1/sec
    K1 = kn1/k1
    K2 = kn2/k2
    K3 = kn3/k3

# Trieber
# Not quite sure why Af, Cf are in units 1/s and order of Mag lower
# since the rxns should still be Ca-dependent
class paramsTrieber:
    k1= 2.16e6 # 1/sec (Af)
    kn1= 400.  # 1/sec
    k2= 30.    # 1/sec (Bf) 
    kn2= 40.   # 1/sec
    k3= 2.57e7 # 1/sec (Cf) 
    kn3= 16.   # 1/sec
    K1 = kn1/k1
    K2 = kn2/k2
    K3 = kn3/k3

class paramsTrieberPLB:
    k1= 2.16e6 # 1/sec (Af)
    kn1= 400.  # 1/sec
    k2= 44.    # 1/sec (Bf) 
    kn2= 50.9e3# 1/sec
    k3= 404e3  # 1/sec (Cf) 
    kn3= 16.   # 1/sec
    K1 = kn1/k1
    K2 = kn2/k2
    K3 = kn3/k3

class paramsSS:            
    K1 = 400/2.16e6 
    K2 = 40/30.         
    K3 = 16/2.57e7

def Inesi(): 
  kn1 = sp.Symbol('kn1')
  k1 = sp.Symbol('k1')
  kn2 = sp.Symbol('kn2')
  k2 = sp.Symbol('k2')
  kn3 = sp.Symbol('kn3')
  k3 = sp.Symbol('k3')
  
  Ca = sp.Symbol('Ca')
  Et = sp.Symbol('Et')
  E = sp.Symbol('E')
  ECa = sp.Symbol('ECa')
  EpCa = sp.Symbol('EpCa')
  EpCa2 = sp.Symbol('EpCa2')
  
  ## slow kinetics 
  # pg 168 of Inesi 1988 paper 
  Eeq = sp.Eq(kn1*kn2*kn3*EpCa2 / (k1*k2*k3*Ca*Ca),E)
  ECaeq = sp.Eq(kn2*kn3*EpCa2 / (k2*k3*Ca),ECa)
  EpCaeq = sp.Eq(kn3*EpCa2 / (k3*Ca),EpCa)
  
  # sub into Et expression
  Eteq = sp.Eq(E + ECa+EpCa + EpCa2,Et)
  Eteqs = Eteq.subs(\
    {E:sp.solve(Eeq,E)[0],\
    ECa:sp.solve(ECaeq,ECa)[0],\
    EpCa:sp.solve(EpCaeq,EpCa)[0]\
    })
  
  sol= sp.solve(Eteqs,EpCa2)[0]
  EpCa2eq = sp.Eq(EpCa2,sol)
  

  EpCa2func = sp.solve(EpCa2eq,EpCa2)[0]
  def func(Cav,p): #p_k1,p_kn1,p_k2,p_kn2,p_k3,p_kn3):
    val    = EpCa2func.subs({\
    Et: 1.,     # assume Et is fixed?? 
    k1: p.k1,   # 1/M sec (Af)
    kn1: p.kn1, # 1/sec
    k2: p.k2,   # 1/sec (Bf) 
    kn2: p.kn2, # 1/sec
    k3: p.k3,   # 1/M sec (Cf) 
    kn3: p.kn3, # 1/sec
    Ca:Cav
    })
    return val.evalf()
  #val  = EpCa2func.subs({Ca:Cav})   

  vfunc = np.vectorize(func,otypes=[np.float])

  return vfunc
  


# functions based on model 
# The K1 terms reflect the equilibrium constants between states, if I recall correctly 
def InesiSS():
  ## stdy kinetics
  K1 = sp.Symbol('K1')
  K2 = sp.Symbol('K2')
  K3 = sp.Symbol('K3')
  Ca = sp.Symbol('Ca')
  Et = sp.Symbol('Et')
  E = sp.Symbol('E')
  ECa = sp.Symbol('ECa')
  EpCa = sp.Symbol('EpCa')
  EpCa2 = sp.Symbol('EpCa2')
  EeqSS = sp.Eq(K1*K2*K3*EpCa2 / (Ca*Ca),E)
  ECaeqSS = sp.Eq(K2*K3*EpCa2 / (Ca),ECa)
  EpCaeqSS = sp.Eq(K3*EpCa2 / (Ca),EpCa)
  
  # sub into Et expression
  Eteq = sp.Eq(E + ECa+EpCa + EpCa2,Et)
  EteqSSs = Eteq.subs(\
    {E:sp.solve(EeqSS,E)[0],\
    ECa:sp.solve(ECaeqSS,ECa)[0],\
    EpCa:sp.solve(EpCaeqSS,EpCa)[0]\
    })
  
  sol= sp.solve(EteqSSs,EpCa2)[0]
  EpCa2SSeq = sp.Eq(EpCa2,sol)

  EpCa2funcSS = sp.solve(EpCa2SSeq,EpCa2)[0]
  print EpCa2funcSS
  def funcSS(Cav,p): #p_k1,p_kn1,p_k2,p_kn2,p_k3,p_kn3):
    val    = EpCa2funcSS.subs({\
    Et: 1.,     # assume Et is fixed?? 
    K1: p.K1,   # 1/M (Af)
    K2: p.K2,   # 1/(Bf) 
    K3: p.K3,   # 1/M(Cf) 
    Ca:Cav
    })
    return val.evalf()
  #val  = EpCa2func.subs({Ca:Cav})   
  vfuncSS = np.vectorize(funcSS,otypes=[np.float])

  return vfuncSS






# Test
def doKin(dfunc,p,cs):
  cEpCa = dfunc(cs,p)   
  #KCa = cs[np.argmin(np.abs(cEpCa-0.5))] * M_to_uM
  KCa = np.interp(0.5,cEpCa,cs)*M_to_uM
  return (cEpCa,KCa)

def doit(fileIn): 

  cs = np.linspace(-8,-4,30)
  cs = 10**cs
  
  vfunc = Inesi()
  vfuncSS = InesiSS()
  
  ## compare models 
  (valsInesi, KCaInesi) = doKin(vfunc,paramsInesi,cs)
  (valsInesiSS, KCaInesiSS) = doKin(vfuncSS,paramsInesi,cs)  
  import matplotlib.pyplot as plt
  plt.figure()
  plt.plot(cs,valsInesi,label="Inesi")
  plt.plot(cs,valsInesiSS,'k.',label="InesiSS") 
  plt.title("SERCA models")         
  plt.xscale('log')
  plt.xlabel('pCa')        
  plt.ylabel('Saturation') 
  plt.legend(loc=0)
  plt.gcf().savefig("sercaModelInesi.png")


# add back P1 state 
# compare w Hig. paper 
import sys
#
# Revisions
#       10.08.10 inception
#

if __name__ == "__main__":
  msg="""
Purpose: 

Usage:
  script.py <arg>

Notes:
"""
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-arg1"):
      arg1=sys.argv[i+1] 




  doit(fileIn)


