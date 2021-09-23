"""Run requests with mpirun"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import subprocess
import sys
import os

#This package
from . import filepath
from . import request
from . import yaml_manager
from . import simultaneous
from . import locators
from . import logging

logger=logging.getLogger(__name__)

#Default tmpfile
DEFAULT_TMPFILE = filepath.Path('.').expanduser().resolve() / 'tmp.yaml'

#Constants from simultaneous
main_mod_name=simultaneous.main_mod_name
work_path=simultaneous.work_path

_MPIRunRequest_props_schema_yaml="""#MPIRunRequest
name: {type: string}
numproc:
  type: integer
  minimum: 1
child: {type: request}
tmpfile: {type: pathlike}
"""

class MPIRunRequest(request.Request):
  """Run the child request using mpirun
  
  User-defined attributes:
  
    - numproc: number of processes to specify with mpirun
    - child: a request to run
    - tmpfile: optional, Path to temporary file to contain the request"""
  _self_task=True #Otherwise, the child won't be run with MPI
  _config_attrs=['numproc','child','tmpfile']
  _validation_schema=request.Request.update_schema(_MPIRunRequest_props_schema_yaml)
  _validation_schema.required=['name','numproc','child']
  _child_attrs=['child']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(MPIRunRequest, self).__init__(**kwargs)
    #tmpfile
    if hasattr(self,'tmpfile'):
      self.tmpfile=filepath.Path(self.tmpfile).expanduser().resolve()
    else:
      self.tmpfile=DEFAULT_TMPFILE
    #Get the input and output files of the child request
    self._more_inputfiles=self.child.inputfiles
    self._more_outputfiles=self.child.outputfiles
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Create the output directories for the request, one time, to avoid conflicts
    self.child.assure_output_dirs()
    #Write the input file
    sdf=locators.SetDataFolder()
    ufs=locators.UpdateFolderStructure()
    yaml_manager.writefile([sdf,ufs,self.child],str(self.tmpfile))
    #Call MPIrun to start the processes
    args=('mpirun','-np','%d'%self.numproc,sys.executable,'-m',main_mod_name,str(self.tmpfile))
    p=subprocess.Popen(args,cwd=work_path,shell=False)
    #Wait for completion
    retcode=p.wait()
    #Clean up
    os.remove(str(self.tmpfile))
    
#Register for loading from yaml
yaml_manager.register_classes([MPIRunRequest])
