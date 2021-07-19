"""Create the input files from templates for run001

Variables that must be set in the mesh template file (.geo):

  - inclusion_dist: distance from pore surface to inclusion

Variables that must be set in the model template file (.yaml):

  - basename
  - meshname
  - modelname
  - sphere_conc_A, sphere_conc_B: species concentrations at the sphere surface
  - top_conc_A, top_conc_B: species concentrations on the "top" reservoir surface
  - bot_conc_A, bot_conc_B: species concentrations on the "bottom" reservoir surface
  - central_pore_pot: electric potential on the surface of the central pore
  - sphere_pot: electric potential on the surface of the sphere
  - debye_length: debye length for the electrostatic potential (1/kappa)

"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from copy import deepcopy
from collections import OrderedDict as odict
from itertools import cycle
import os
import os.path as osp
import sys

#Site packages
import yaml
from jinja2 import Template

#Locate source folders
if 'SRCLOC' in os.environ.keys():
  srcfolder=osp.normpath(osp.abspath(os.environ['SRCLOC']))
else:
  srcfolder=osp.abspath(osp.split(__file__)[0])

#add python code folder(s) to path
if not srcfolder in sys.path:
  sys.path.append(srcfolder)

import folderstructure as FS

basename='run003'

#Constant model parameters
model_consts={'sphere_conc_A': 0.0,
              'sphere_conc_B': 'null',
              'top_conc_A': 6.0e-4,
              'top_conc_B': 0.0,
              'bot_conc_A': 0.0,
              'bot_conc_B': 0.0,
              'sphere_pot': -25.0e-3}

geo_template_file = osp.join(FS.datafolder,"mesh/templates/%s.geo.jinja2"%basename)
mesh_yaml_outfile=osp.join(FS.params_mesh_folder,'%s.yaml'%basename)
geo_outdir=osp.join(FS.geofolder,basename)
meta_outdir=osp.join(FS.meshmeta_outfolder,basename)
for dirname in [geo_outdir, meta_outdir]:
  if not osp.isdir(dirname):
    os.mkdir(dirname)
model_template_file = osp.join(FS.datafolder,"model/templates/model_%s.yaml.jinja2"%basename)
model_yaml_outfile_tmpl=osp.join(FS.params_model_folder,basename+'_%s.yaml')
simulation_script_out=basename+'_sims.sh'
codebook_outfile='codebook_%s.yaml'%basename

dist_list = [0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
pore_pot_list = [0.0]
kappa_list=[1.0]
##dist_list = [0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
##pore_pot_list = [0.0, 25.0e-3,-25.0e-3]
##kappa_list=[0.01, 0.1, 1.0, 10.0, 100.0]

#The number of characters in the list below is the number of runs done in parallel
runlist='ABCDE'


def load_template(tmpl_fpath):
  """Load the jinja2 template from the given file.
  
  Note that this method does not support extending templates."""
  with open(tmpl_fpath,'r') as fh:
    tdata=fh.read()
  return Template(tdata)

def apply_template(template,out_fpath,values):
  """Apply the given data to the specified template to produce an output file.
  
  Arguments:
  
    - template = jinja2 Template instance
    - out_fpath = path to output file, as string
    - values = dictionary of values to put into the template, {template variable name: value}"""
  #Substitute output
  out_data=template.render(**values)
  #Write output file
  with open(out_fpath,'w') as fh:
    fh.write(out_data)

#Read the yaml input templates
geo_tmpl=load_template(geo_template_file)
model_tmpl=load_template(model_template_file)
with open(geo_template_file,'r') as fh:
  geo_tmpl_data=fh.read()
with open(model_template_file,'r') as fh:
  model_tmpl_data=fh.read()

#Model yaml outfiles for doing simultaneous runs
model_yaml_outfile_list=[]
for runcode in runlist:
  model_yaml_outfile_list.append(model_yaml_outfile_tmpl%runcode)
model_yaml_outfiles=cycle(model_yaml_outfile_list)

#Initialize yaml output
output_files=odict()
output_files[codebook_outfile]="""%YAML 1.2\n---\n"""
output_files[mesh_yaml_outfile]="""%YAML 1.2\n"""
for myo in model_yaml_outfile_list:
  output_files[myo]="""%YAML 1.2\n"""
output_files[simulation_script_out]="#!/bin/bash\n#Run simultaneous simulations\n\n"

#Script for simultaneous runs
for runcode in runlist:
  output_files[simulation_script_out]+='python2 $SRCLOC/simulator_run.py params/model/%s_%s.yaml >logs/%s_%s.txt 2>&1 &\n'%(basename,runcode,basename,runcode)

modelid=0
#Loop over meshes
#Only one varying parameter: inclusion distance
for meshid,dist in enumerate(dist_list):
  #Create the file from the template
  mesh_values={'inclusion_dist':dist}
  meshname='%s_mesh_%03d'%(basename,meshid)
  geo_outfpath=osp.join(geo_outdir,'%s.geo'%meshname)
  apply_template(geo_tmpl,geo_outfpath,mesh_values)
  #Create a basic mesh metadata file
  meta_str='\n'.join(['%s: %s'%tup for tup in mesh_values.items()])
  with open(osp.join(meta_outdir,'%s.yaml'%meshname),'w') as fh:
    fh.write(meta_str)
  #Add to the yaml output
  yaml_doc="---\nmeshname: %s\ngeomdefname: null\ntmplvalues:\n"%(meshname)
  yaml_doc+="\n".join(["  %s: %s"%tup for tup in mesh_values.items()])+"\n"
  output_files[mesh_yaml_outfile]+=yaml_doc
  #Loop over simulations to be done with this mesh
  for pore_pot in pore_pot_list: #Each pore potential
    #A potential of zero means no dependence on kappa
    for kappa in kappa_list: #Each kappa value
      modelname='%s_model_%03d'%(basename,modelid)
      #Add this to the codebook
      output_files[codebook_outfile]+="%s:\n  inc_dist: %s\n  pore_pot: %s\n  kappa: %s\n"%(modelname,dist,pore_pot,kappa)
      #Get the dictionary of values for the template
      model_values=deepcopy(model_consts)
      model_values['basename']=basename
      model_values['meshname']=meshname
      model_values['modelname']=modelname
      model_values['central_pore_pot']=pore_pot
      model_values['debye_length']=1.0/kappa
      #Add model to yaml document
      yaml_doc=model_tmpl.render(**model_values)
      output_files[next(model_yaml_outfiles)]+=yaml_doc
      #Next model
      modelid+=1


#Write the output yaml files
for fpath,data in output_files.items():
  with open(fpath,"w") as fh:
    fh.write(data)

