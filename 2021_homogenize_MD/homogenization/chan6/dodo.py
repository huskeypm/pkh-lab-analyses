"""Doit file for running requests

Basic Usage:
doit
For more information, try
doit help

To specify an alternate control file,
doit control=controlfile"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#Local
from simproc import *
from simproc.requesthandler.cmdline import yield_doit_tasks

def task_all():
  return yield_doit_tasks()
