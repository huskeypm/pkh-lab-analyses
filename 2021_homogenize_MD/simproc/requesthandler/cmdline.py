"""Command-line support for running requests, and returning requests as doit tasks"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import doctest
import importlib
import sys

#Site packages
from doit import get_var

#Local
from . import *

#Constants
default_controlfile = locators.DATAFOLDER / 'control.yaml'

def run_validation(verbose=False):
  """Run the validation tests"""
  reslist=[]
  for m in doctest_modules:
    print("Testing module %s."%m.__name__)
    reslist.append(doctest.testmod(m,verbose=verbose))
    print("---")
  for fpath in doctest_files:
    print("Running tests in %s."%fpath)
    reslist.append(doctest.testfile(fpath,module_relative=False,verbose=verbose)) #To make sure the file can be found when the package is zipped, we have already found its absolute path
    print("---")
  fails,atts=[sum(l) for l in zip(*reslist)]
  print("Passed %d/%d total"%(atts-fails,atts))

def run():
  """Run request files from the command line"""
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('requestfile',nargs="*",help="Path to file containing the request(s) to run. Multiple request files may be listed.")
  parser.add_argument('--verbose',action='store_true',help="Provide verbose output where appropriate.")
  parser.add_argument('--tasks_only',action='store_true',help="Run only requests that define tasks.")
  parser.add_argument('--modules',nargs="+",metavar="MODULE",help="Additional python modules defining classes loadable from yaml input")
  parser.add_argument('--validate',action='store_true',help="Perform validation. If requestfiles are also listed, validation is run first.")
  parser.add_argument('--select',nargs="+",metavar="TASKNAME",help="Specify which tasks to run, by name.")
  parser.add_argument('--stdout_level',nargs="?",help="Minimum level for events to be logged to stdout. If not provided, no events logged to stdout.",default=None)
  parser.add_argument('--stderr_level',nargs="?",help="Minimum level for events to be logged to stderr. If not provided, no events logged to stderr.",default=None)
  parser.add_argument('--loglevel',nargs="?",help="Minimum level for events to be sent to the log file. If not provided, no log file will be created.",default=None)
  parser.add_argument('--logstem',nargs="?",help="Stem name of the log file.",default="simproc")
  parser.add_argument('--logdir',nargs="?",help="Path to the folder containing the log file. Overrides logdir_rel if given.",default=None)
  parser.add_argument('--logdir_rel',nargs="?",help="Path to the folder containing the log file, relative to the DATAFOLDER.",default="logs")
  parser.add_argument('--logext',nargs="?",help="Extension for the log file.",default=".log.yaml")
  parser.add_argument('--log_num_digits',nargs="?",help="Number of digits in the ID portion of the log file name.",default=3)
  parser.add_argument('--log_sepchar',nargs="?",help="Separation character within the log filename before the ID.",default=".")
  cmdline=parser.parse_args()
  
  #run validation if requested
  if cmdline.validate:
    run_validation(verbose=cmdline.verbose)

  #Set up logging
  if cmdline.loglevel is not None:
    logging.configure_logfile(level=cmdline.loglevel,stem=cmdline.logstem,logdir_rel=cmdline.logdir_rel,
                              logdir_abs=cmdline.logdir,ext=cmdline.logext,
                              num_digits=cmdline.log_num_digits,sepchar=cmdline.log_sepchar)
  if cmdline.stdout_level is not None:
    logging.configure_stdout(cmdline.stdout_level)
  if cmdline.stderr_level is not None:
    logging.configure_stderr(cmdline.stderr_level)

  #Confirm that specified request file(s) exist(s)
  file_list=[filepath.Path(rf,isFile=True) for rf in cmdline.requestfile]
  nonexist=[rf for rf in file_list if not rf.exists()]
  assert len(nonexist)==0, "Could not find specified request file(s) %s"%str([rf.fullpath for rf in nonexist])

  #Load the requested modules
  if cmdline.modules is not None:
    customization.load_modules(cmdline.modules)

  #Initialize a RequestFileListRequest
  req=requestfile.RequestFileListRequest(name="(from command line)",requestfiles=file_list)

  #Run
  if cmdline.tasks_only or cmdline.select is not None:
    for treq in req.all_task_requests():
      if cmdline.select is None or getattr(treq,'name','') in cmdline.select:
        if cmdline.verbose:
          print(treq.name)
        treq.run()
  else:
    req.run()

def yield_doit_tasks():
  """Create task definitions for doit
  
  Usage:
    ```def task_all():
        return yield_doit_tasks()```
  """
  #Get the controlfile as a Path
  controlfile=get_var('control',default_controlfile)
  controlfile=filepath.Path(controlfile)

  #Initialize a RequestFileRequest
  req=requestfile.RequestFileRequest(requestfile=controlfile)

  #Return the task generator
  return req.all_tasks()

