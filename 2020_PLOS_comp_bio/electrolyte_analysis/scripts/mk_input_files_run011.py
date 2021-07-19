"""Create the input files from templates for run009

Variables that must be set in the mesh template file (.geo):

  - pore_radius = radius of pore
  - RevW = width of reservoir (X)
  - RevH = height of reservoir (Y)
  - RevD = depth of reservoir (Z)
  - pore_length = length of pore
  - inclusionlist = = sequence of inclusions, each inclusion a tuple (id, x,y,z,r)
  - mcar1 = coarsest scale, for ends of reservoirs
  - mcar2 = for porous surface away from pore
  - mcar3 = for pore boundaries
  - mcar4 = finest scale, for spherical inclusion

Variables that must be set in the model template file (.yaml):

  - basename
  - meshname
  - modelname
  - z_A, z_B, z_C: species charges
  - top_conc_A, top_conc_B, top_conc_C: species concentrations on the upper reservoir surface
  - bot_conc_A, bot_conc_B, bot_conc_C: species concentrations on the lower reservoir surface
  - pore_conc_A, pore_conc_B, pore_conc_C: species concentrations on the pore surface
  - top_potential, bot_potential: electric potential values at the upper and lower reservoir face
  - debye_length: debye length for the electrostatic potential (1/kappa)
  - r_inc_conc_A, r_inc_conc_B: concentration values for Dirichlet boundary conditions on the reactive inclusions
  - nr_inc_conc_A, nr_inc_conc_B: concentration values for Dirichlet boundary conditions on the nonreactive inclusions
  - r_inc_potential: electric potential for Dirichlet boundary condition on the reactive inclusions
  - nr_inc_potential: electric potential for Dirichlet boundary condition on the nonreactive inclusions
  - all_inc_ids: list of all inclusion id values
  - reactiveAB_inc_ids: list of id values for only reactive A->B inclusions
  - reactiveBC_inc_ids: list of id values for only reactive B->C inclusions
  - nonreactive_inc_ids: list of id values for only non-reactive inclusions

Variables controlling inclusion generation:

  - randseed = random seed value for inclusion generation, use -1 for regular lattice

"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from argparse import Namespace
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

basename='run011'

#The number of characters in this string is the number of runs done in parallel
runlist='A'

#Inclusions (must be in top half of pore): (ID,x,y,z,rad)
inclusionlist=[(1000,0.0,0.0,35.0,1.25),(1100,0.0,0.0,-30.0,1.25)]

#Parameter values: each combination of entries is used
mesh_vars=odict([
  ('pore_radius', [3.0]),
  ('RevW', [12.0]),
  ('RevH', [12.0]),
  ('RevD', [40.0]),
  ('pore_length', [90.0]),
  ('inclusionlist',[inclusionlist]),
  ('mscale', [1.0]),
  ('mcar1', [5.0]),
  ('mcar2', [2.0]),
  ('mcar3', [1.0]),
  ('mcar4', [0.2])]
)

model_vars=odict([
  ('z_A', [1]),
  ('z_B', [1]),
  ('z_C', [1]),
  ('num_reactiveAB_inclusions', [1]),
  ('num_reactiveBC_inclusions', [1]),
  ('top_conc_A', [6.0e-4]),
  ('top_conc_B', [0.0]),
  ('top_conc_C', [0.0]),
  ('bot_conc_A', [0.0]),
  ('bot_conc_B', [0.0]),
  ('bot_conc_C', [0.0]),
  ('pore_conc_A', ['null']),
  ('pore_conc_B', ['null']),
  ('pore_conc_C', ['null']),
  ('top_potential', [0.0]),
  ('bot_potential', [0.0]),
  ('rAB_inc_conc_A', [0.0]),
  ('rAB_inc_conc_B', [1.0e-4]),
  ('rBC_inc_conc_B', [0.0]),
  ('rBC_inc_conc_C', [1.0e-4]),
  ('rAB_inc_potential', [-25.0e-3]),
  ('rBC_inc_potential', [-25.0e-3]),
  ('nr_inc_conc_A', ['null']),
  ('nr_inc_conc_B', ['null']),
  ('nr_inc_conc_C', ['null']),
  ('nr_inc_potential', [0]),
  ('central_pore_potential', [0.0]),
  ('kappa', [1.0])]
)

#  ('', []),

mesh_direct_fields=['pore_radius','pore_length','RevW','RevH','RevD',
                    'mscale','mcar1','mcar2','mcar3','mcar4','inclusionlist']

model_direct_fields=['basename','meshname','modelname',
                    'z_A','z_B','z_C',
                    'top_conc_A','top_conc_B','top_conc_C',
                    'bot_conc_A','bot_conc_B','bot_conc_C',
                    'pore_conc_A','pore_conc_B','pore_conc_C',
                    'top_potential','bot_potential',
                    'rAB_inc_conc_A','rAB_inc_conc_B',
                    'rBC_inc_conc_B','rBC_inc_conc_C',
                    'rAB_inc_potential','rBC_inc_potential',
                    'nr_inc_conc_A','nr_inc_conc_B','nr_inc_conc_C','nr_inc_potential','central_pore_potential']

model_direct_mesh_fields=['all_inc_ids']  

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
output_exclusions=['left_inclusions','right_inclusions','left_ids','right_ids']

def load_template(tmpl_fpath):
  """Load the jinja2 template from the given file.
  
  Note that this method does not support extending templates."""
  with open(tmpl_fpath,'r') as fh:
    tdata=fh.read()
  return Template(tdata,trim_blocks=True)

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

def output_prep(d,exclusions=None):
  """Convert dictionary to sequence of tuples, with proper string escaping, and excluding select keys"""
  if exclusions is None:
    exclusions = []
  outlist=[]
  for k,v in d.items():
    if k not in exclusions:
      if isinstance(v,str):
        out='"%s"'%str(v)
      else:
        out="%s"%str(v)
      outlist.append((k,out))
  return outlist

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

def calc_tmpl_mesh(mesh_values):
  """Compute the values that go into the mesh template, from the relevant inputs"""
  output={}
  #Direct copy
  for k in mesh_direct_fields:
    output[k]=mesh_values[k]
  output['top_inclusion_ids']="".join([',%d'%t[0] for t in mesh_values['inclusionlist'] if t[3]>0.0])
  output['bot_inclusion_ids']="".join([',%d'%t[0] for t in mesh_values['inclusionlist'] if t[3]<0.0])
  output['all_inc_ids']=[t[0] for t in mesh_values['inclusionlist']]
  return output

def calc_tmpl_model(mesh_data, model_values):
  """Compute the values that go into the model template, from the relevant inputs"""
  output={}
  #Direct copy
  for k in model_direct_fields:
    output[k]=model_values[k]
  for k in model_direct_mesh_fields:
    output[k]=mesh_data[k]
  #Calculations
  output['debye_length']=1.0/model_values['kappa']
  #Separate reactive and non-reactive inclusions
  nrAB=model_values['num_reactiveAB_inclusions']
  nrBC=model_values['num_reactiveBC_inclusions']
  nr=nrAB+nrBC
  output['reactiveAB_inc_ids']=mesh_data['all_inc_ids'][:nrAB]
  output['reactiveBC_inc_ids']=mesh_data['all_inc_ids'][nrAB:nr]
  output['nonreactive_inc_ids']=mesh_data['all_inc_ids'][nr:]
  return output

#Loop over meshes
modelid = 0
meshvals_iter=product(*list(mesh_vars.values()))
for meshid,meshvals in enumerate(meshvals_iter):
  meshname='%s_mesh_%03d'%(basename,meshid)
  mesh_values=odict(zip(mesh_vars.keys(),meshvals))
  #Compute template input
  mesh_tmpl_data=calc_tmpl_mesh(mesh_values)
  #Create the file from the template
  geo_outfpath=osp.join(geo_outdir,'%s.geo'%meshname)
  apply_template(geo_tmpl,geo_outfpath,mesh_tmpl_data)
  #Add to the yaml output
  yaml_doc="---\nmeshname: %s\ngeomdefname: null\ntmplvalues:\n"%(meshname)
  yaml_doc+="\n".join(["  %s: %s"%tup for tup in output_prep(mesh_values,output_exclusions)])+"\n"
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
    #Compute template input
    model_tmpl_data=calc_tmpl_model(mesh_tmpl_data,model_values)
    #Add model to yaml document
    yaml_doc=model_tmpl.render(**model_tmpl_data)
    output_files[next(model_yaml_outfiles)]+=yaml_doc
    #Add model to the codebook
    all_values=odict()
    all_values.update(mesh_values)
    all_values.update(model_values)
    code_entry="%s:\n"%modelname
    for pair in output_prep(all_values,output_exclusions):
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

