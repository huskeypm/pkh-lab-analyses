"""Interface for dealing with pickle files"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import pickle

#Constants
##pickle_protocol = 4 #The newest protocol, requires python 3.4 or above.
pickle_protocol = 2 #For compatibility with Python 2

def readfile(fpath):
  "read object from pickle file"
  with open(str(fpath), 'rb') as fp:
    obj=pickle.load(fp)
  return obj

def writefile(obj,fpath):
  "write object to pickle file, overwriting"
  with open(str(fpath), 'wb') as fp:
    pickle.dump(obj,fp,pickle_protocol)
  return
