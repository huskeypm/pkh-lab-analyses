"""Define the expected folder structure"""

#Standard library
import os
import os.path as osp
import sys

#Locate source folder
if 'SRCLOC' in os.environ.keys():
  srcfolder=osp.normpath(osp.abspath(os.environ['SRCLOC']))
else:
  srcfolder=osp.abspath(osp.split(__file__)[0])

#Add the source folder to the python path if not already present,
#to allow importing modules
#add python code folder(s) to path
if not srcfolder in sys.path:
  sys.path.append(srcfolder)

from filepath import Path

#Use filepath.Path now
srcfolder=Path(srcfolder,isFile=False)

# simulator_modules_folder=srcfolder / 'simulators'
# if not simulator_modules_folder in sys.path:
#   sys.path.append(simulator_modules_folder)

#Locate data folder
if 'DATALOC' in os.environ.keys():
  datafolder=Path(osp.normpath(osp.abspath(os.environ['DATALOC'])))
else:
  datafolder=srcfolder.parent / 'data'
