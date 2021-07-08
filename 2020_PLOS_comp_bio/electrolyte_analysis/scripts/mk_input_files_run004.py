"""Create the input files from templates for run004

Variables that must be set in the mesh template file (.geo):

  - pore_radius = radius of pore
  - RevW = width of reservoir (X)
  - RevH = height of reservoir (Y)
  - RevD = depth of reservoir (Z)
  - pore_length = length of pore
  - incRadius = radius of spherical inclusion
  - incDist = distance from pore surface to inclusion
  - incDepth =  how far inside the pore the inclusion is
  - mcar1 = coarsest scale, for ends of reservoirs
  - mcar2 = for porous surface away from pore
  - mcar3 = for pore boundaries
  - mcar4 = finest scale, for spherical inclusion

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
from itertools import cycle, product
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

basename='run004'

#The number of characters in this string is the number of runs done in parallel
runlist='ABCDE'

#Parameter values: each combination of entries is used
mesh_vars=odict([
  ('pore_radius', [3.0]),
  ('RevW', [12.0]),
  ('RevH', [12.0]),
  ('RevD', [40.0]),
  ('pore_length', [90.0]),
  ('incRadius', [1.5]),
  # ('incDist', [0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
  ('incDist', [0.5]),
  ('incDepth',  [2.0]),
  ('mscale', [1.0]),
  ('mcar1', [5.0]),
  ('mcar2', [2.0]),
  ('mcar3', [1.0]),
  ('mcar4', [0.2])]
)

model_vars=odict([
  ('sphere_conc_A', [0.0]),
  ('sphere_conc_B', ['null']),
  ('top_conc_A', [6.0e-4]),
  ('top_conc_B', [0.0]),
  ('bot_conc_A', [0.0]),
  ('bot_conc_B', [0.0]),
  ('central_pore_pot', [0.0]),
  ('sphere_pot', [-25.0e-3]),
  ('kappa', [1.0])]
)

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
myo_start="""%YAML 1.2\n"""
output_files=odict()
output_files[codebook_outfile]="""%YAML 1.2\n---\n"""
output_files[mesh_yaml_outfile]="""%YAML 1.2\n"""
for myo in model_yaml_outfile_list:
  output_files[myo]=myo_start
output_files[simulation_script_out]="#!/bin/bash\n#Run simultaneous simulations\n\n"

#Loop over meshes
modelid = 0
meshvals_iter=product(*list(mesh_vars.values()))
for meshid,meshvals in enumerate(meshvals_iter):
  #Create the file from the template
  meshname='%s_mesh_%03d'%(basename,meshid)
  mesh_values=odict(zip(mesh_vars.keys(),meshvals))
  geo_outfpath=osp.join(geo_outdir,'%s.geo'%meshname)
  apply_template(geo_tmpl,geo_outfpath,mesh_values)
  #Add to the yaml output
  yaml_doc="---\nmeshname: %s\ngeomdefname: null\ntmplvalues:\n"%(meshname)
  yaml_doc+="\n".join(["  %s: %s"%tup for tup in mesh_values.items()])+"\n"
  output_files[mesh_yaml_outfile]+=yaml_doc
  #Loop over simulations to be done with this mesh
  modelvals_iter=product(*list(model_vars.values()))
  for modelvals in modelvals_iter:
    modelname='%s_model_%03d'%(basename,modelid)
    model_values=odict()
    model_values['basename']=basename
    model_values['meshname']=meshname
    model_values['modelname']=modelname
    model_values.update(odict(zip(model_vars.keys(),modelvals)))
    if model_values['kappa']==0:
      model_values['debye_length']='null'
    else:
      model_values['debye_length']=1.0/model_values['kappa']
    #Add model to yaml document
    yaml_doc=model_tmpl.render(**model_values)
    output_files[next(model_yaml_outfiles)]+=yaml_doc
    #Add model to the codebook
    all_values=odict()
    all_values.update(mesh_values)
    all_values.update(model_values)
    code_entry="%s:\n"%modelname
    for pair in all_values.items():
      code_entry+="  %s: %s\n"%pair
    output_files[codebook_outfile]+=code_entry
    #Next model
    modelid += 1

#Eliminate model yaml output files with no data
final_model_yaml_outfile_list=[]
for myo in model_yaml_outfile_list:
  if len(output_files[myo])>len(myo_start):
    final_model_yaml_outfile_list.append(myo)
  else:
    junk=output_files.pop(myo)

#Script for simultaneous runs
for runcode in runlist:
  if model_yaml_outfile_tmpl%runcode in final_model_yaml_outfile_list:
    output_files[simulation_script_out]+='python2 $SRCLOC/simulator_run.py params/model/%s_%s.yaml >logs/%s_%s.txt 2>&1 &\n'%(basename,runcode,basename,runcode)

#Write the output yaml files
for fpath,data in output_files.items():
  with open(fpath,"w") as fh:
    fh.write(data)

