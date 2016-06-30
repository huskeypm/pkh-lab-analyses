
import homoglight as hl 
import numpy as np


def doit():                 
  fileXML="example/volfracs/volFrac_0.50_mesh.xml.gz"
  problem = hl.runHomog(fileXML,verbose=True)
  assert(np.abs(problem.d_eff[0]-0.3939)<0.01), "Don't commit! somthing changed"
  #assert(np.abs(4-0.3939)<0.01), "Don't commit! somthing changed"
  print "All is ok!"
  quit()

#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Very simple consistency check 
 
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
      doit()      
  





  raise RuntimeError("Arguments not understood")




