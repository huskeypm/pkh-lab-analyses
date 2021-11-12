"""Allow command-line execution"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#This Package
from . import *
from .requesthandler import cmdline

#Handle command-line execution
if __name__ == '__main__':
  cmdline.run()

