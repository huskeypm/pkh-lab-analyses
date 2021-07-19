"""Define the expected folder structure"""

#Standard library
import os
import os.path as osp
import sys

#Locate source folders
if 'SRCLOC' in os.environ.keys():
  srcfolder=osp.normpath(osp.abspath(os.environ['SRCLOC']))
else:
  srcfolder=osp.abspath(osp.split(__file__)[0])

#subfolders of src
custom_modules_folder=osp.join(srcfolder,'customizations')
simulator_modules_folder=osp.join(srcfolder,'simulators')

#add python code folder(s) to path
if not srcfolder in sys.path:
  sys.path.append(srcfolder)
if not custom_modules_folder in sys.path:
  sys.path.append(custom_modules_folder)
if not simulator_modules_folder in sys.path:
  sys.path.append(simulator_modules_folder)

#Locate data folder
if 'DATALOC' in os.environ.keys():
  datafolder=osp.normpath(osp.abspath(os.environ['DATALOC']))
else:
  datafolder=osp.join(osp.split(srcfolder)[0],'data')

#params
paramsfolder=osp.join(datafolder,'params')
params_mesh_folder=osp.join(paramsfolder,'mesh')
params_model_folder=osp.join(paramsfolder,'model')
params_postproc_folder=osp.join(paramsfolder,'postproc')
params_paramgen_folder=osp.join(paramsfolder,'paramgen')

#mesh
meshfolder=osp.join(datafolder,'mesh')
geomdef_folder=osp.join(meshfolder,'geomdef')
geotemplates_folder=osp.join(meshfolder,'templates')
geofolder=osp.join(meshfolder,'geo')
mshfolder=osp.join(meshfolder,'msh')
xmlfolder=osp.join(meshfolder,'xml')
mesh_hdf5_folder=osp.join(meshfolder,'hdf5')
gmsh_outfolder=osp.join(meshfolder,'gmsh_out')
meshmeta_outfolder=osp.join(meshfolder,'metadata')

#solutions
solnfolder=osp.join(datafolder,'solutions')
infofile="info.yaml" #Output file describing the run, including select results

#post-processing
postprocfolder=osp.join(datafolder,'postproc')

#parameter generation
pgtemplates_folder=osp.join(datafolder,'paramgen_tmpl')
