"""For running simulations from the input yaml file(s)"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import os.path as osp
import sys

#Site packages

#Local
import folderstructure as FS
import common
import buildgeom
import simulator_general

#Simulator modules
simulator_module_list=['projector','fickian_unhomog','smol_unhomog',
                       'smol_reactive_surface','smol_reactive_surface_subdomains','smol_reactive_surface_firstorder',
                       'tdpnp_unhomog','fickian_homog']

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

class ModelParameters(simulator_general.ModelParametersBase):
  """Extend ModelParametersBase to allow for task generation and execution of arbitrary simulator module

  Attributes:

    - simulatorclass = the class used to run the simulation"""
  __slots__=('simulatorclass',)
  #Load the simulator modules and map equation names to their simulator classes
  #Mapping from ModelParameters.equation to the appropriate simulator classes
  #Each module implementing a simulator should define a simulatorclasses dictionary, so we just need to put them all together
  simulatorclasses={}
  for sm_name in simulator_module_list:
    simulatorclasses.update(importlib.import_module(sm_name).simulatorclasses)
  
  def __init__(self,**kwd):
    #Initialization from base class
    super(ModelParameters, self).__init__(**kwd)
    #Get the simulator class
    self.simulatorclass=self.simulatorclasses[self.equation]
    #Add to code file dependencies
    self._more_inputfiles+=[thisfile, sys.modules[self.simulatorclass.__module__].__file__]

  def run(self):
    """Run the loaded simulation."""
    #Run the simulation
    print(self.modelname)
    self.simulatorclass.complete(self)

#Support command-line arguments
if __name__ == '__main__':
  program_description='Solve a diffusion equation with fenics'
  input_file_description='Path to file containing ModelParameters definitions'
  other_selection={'equation':ModelParameters.simulatorclasses.keys()}
  
  common.run_cmd_line(program_description,input_file_description,ModelParameters,other_selection=other_selection)
  #other_selection is needed so we only try to run models whose equation we have a simulator for
