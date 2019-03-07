#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#import analyzeODE as ao
import analyzeGotran as anG

#
# ROUTINE  
#

# rate - downsample rate
def downsample(fileName, fileOutName=None,rate=10): 
  dataFull = anG.readPickle(fileName)
   
  p = dataFull['p']
  p_idx = dataFull['p_idx']
  s = dataFull['s']; #sDs = s[::rate,]
  s_idx = dataFull['s_idx']
  j = dataFull['j']; #jDs = j[::rate,]
  j_idx = dataFull['j_idx']
  t = dataFull['t']; #tDs = t[::rate]

  sDs,jDs,tDs = downsampleData(s,j,t,rate)
  
  if fileOutName==None:
    redFileName = fileName.replace(".pickle","_red.pickle")
  else: 
    redFileName = fileOutName 
  anG.writePickle(redFileName,p,p_idx,sDs,s_idx,jDs,j_idx,tDs)

import numpy as np 
def downsampleData(s,j,t,rate):
  if isinstance(rate,float):
    print "WARNING: changing rate %f into int"
    rate = np.int( rate ) 
  sDs = s[::rate,]
  jDs = j[::rate,]
  tDs = t[::rate,]
  return sDs, jDs, tDs



#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Downsample pickle file 
 
Usage:
"""
  msg+="  %s -pickleName NAME <-rate INT>" % (scriptName)
  msg+="""
  
 
Notes:
  High-frequency information (below downsample limit) will be lost!!

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
  fileName = False 
  fileOutName = None
  rate = 10 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-pickleName"):
      fileName=sys.argv[i+1] 
    if(arg=="-pickleOutName"):
      fileOutName=sys.argv[i+1] 
    if(arg=="-rate"):
      rate =int(sys.argv[i+1])
  





  if fileName==False:
    raise RuntimeError("Arguments not understood")
  else:
    downsample(fileName,fileOutName=fileOutName, rate=rate)




