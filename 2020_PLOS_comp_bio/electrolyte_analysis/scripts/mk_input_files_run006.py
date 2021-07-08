"""Create the input files from templates for run006

Variables that must be set in the mesh template file (.geo):

  - S = spacing between slots
  - d = width of slot
  - H = reservoir length (distance from slot to opposite end of reservoir)
  - slot_length = length of the slot
  - mcar1 = coarsest scale, for ends of reservoirs
  - mcar2 = for porous surface away from pore
  - mcar3 = for pore boundaries
  - mcar4 = finest scale, for inclusion surfaces
  - left_inclusions, right_inclusions = sequences of inclusions, each inclusion a tuple (id, x,y,r) (calculated, see below)
  - left_ids, right_ids = strings with inclusion ids for defining the surfaces (calculated, see below)

Variables that must be set in the model template file (.yaml):

  - basename
  - meshname
  - modelname
  - left_conc_A, left_conc_B: species concentrations on the left reservoir surface
  - right_conc_A, right_conc_B: species concentrations on the right reservoir surface
  - left_potential, right_potential: electric potential values at the left and right reservoir face
  - debye_length: debye length for the electrostatic potential (1/kappa)
  - r_inc_conc_A, r_inc_conc_B: concentration values for Dirichlet boundary conditions on the reactive inclusions
  - nr_inc_conc_A, nr_inc_conc_B: concentration values for Dirichlet boundary conditions on the nonreactive inclusions
  - r_inc_potential: electric potential for Dirichlet boundary condition on the reactive inclusions
  - nr_inc_potential: electric potential for Dirichlet boundary condition on the nonreactive inclusions
  - all_inc_ids: list of all inclusion id values
  - reactive_inc_ids: list of id values for only reactive inclusions
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
import random_inclusions_2D

basename='run006'

#The number of characters in this string is the number of runs done in parallel
runlist='ABCDEFGHI'

#Parameter values: each combination of entries is used
mesh_vars=odict([
  # ('randseed', [0, 1000]),
  ('randseed', [1005]),
  ('S', [12.0]),
  # ('d', [4.0, 5.0, 6.0]),
  ('d', [5.0]),
  ('H', [40.0]),
  ('slot_length', [90.0]),
  ('n_both', [0,1,2,5,10,25,50,60,70,80,90,100]),
  # ('n_left', [20]),
  # ('n_right', [21]),
  ('mcar1', [5.0]),
  ('mcar2', [2.0]),
  ('mcar3', [1.0]),
  ('mcar4', [0.2])]
)

model_vars=odict([
  ('left_conc_A', [6.0e-4]),
  ('left_conc_B', [0.0]),
  ('right_conc_A', [0.0]),
  ('right_conc_B', [0.0]),
  ('left_potential', [0.0]),
  ('right_potential', [0.0]),
  ('r_inc_conc_A', [0.0]),
  ('r_inc_conc_B', ['null']),
  ('r_inc_potential', [-25.0e-3]),
  ('nr_inc_conc_A', ['null']),
  ('nr_inc_conc_B', ['null']),
  ('nr_inc_potential', [0]),
  ('kappa', [1.0])]
)

#  ('', []),

inc_id_start=1000
maxtries=1e4
class InclusionData(object):
  def format_inclusions(self,centers,rvals):
    global inc_id
    inclusions=[]
    for i,cen in enumerate(centers):
      inclusions.append((inc_id,cen[0],cen[1],rvals[i]))
      inc_id+=10
    return inclusions
  def gen_periodic_inclusions(self):
    bbox=self.box
    centers,rvals=random_inclusions_2D.nonrandom_circles(self.fixed,self.ncirc,
                                                      self.radius_limits[1],
                                                      bbox,self.eps)
    return self.format_inclusions(centers,rvals)
  def gen_random_inclusions(self):
    bbox=self.box
    centers,rvals=random_inclusions_2D.random_circles(self.fixed,self.ncirc,
                                                      (bbox.xmin,bbox.xmax),(bbox.ymin,bbox.ymax),self.radius_limits,
                                                      bbox,maxtries,self.eps)
    return self.format_inclusions(centers,rvals)

class LeftInclusions(InclusionData):
  def __init__(self,slot_length,d,ncirc):
    self.box=random_inclusions_2D.Rect(-slot_length/2,-d/2,0,d/2)
    self.fixed=[(-slot_length*3/8,0,1)]
    self.ncirc=ncirc
    self.radius_limits=(0.25,0.5)
    self.eps=0.25

class RightInclusions(InclusionData):
  def __init__(self,slot_length,d,ncirc):
    self.box=random_inclusions_2D.Rect(0,-d/2,slot_length/2,d/2)
    self.fixed=[]
    self.ncirc=ncirc
    self.radius_limits=(0.25,0.5)
    self.eps=0.25
    

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
  for k in ['S','d','H','slot_length','mcar1','mcar2','mcar3','mcar4']:
    output[k]=mesh_values[k]
  #Number of inclusions on each side
  if 'n_left' not in mesh_values.keys() or 'n_right' not in mesh_values.keys():
    mesh_values['n_left']=mesh_values['n_both']
    mesh_values['n_right']=mesh_values['n_both']
  #Generate inclusions
  global inc_id
  inc_id=inc_id_start
  incdata_l=LeftInclusions(mesh_values['slot_length'],mesh_values['d'],mesh_values['n_left'])
  incdata_r=RightInclusions(mesh_values['slot_length'],mesh_values['d'],mesh_values['n_right'])
  if mesh_values['randseed'] == -1:
    #Use regular lattice of inclusions
    left_inclusions=incdata_l.gen_periodic_inclusions()
    right_inclusions=incdata_r.gen_periodic_inclusions()
  else:
    #Random
    random_inclusions_2D.reseed(mesh_values['randseed'])
    left_inclusions=incdata_l.gen_random_inclusions()
    right_inclusions=incdata_r.gen_random_inclusions()
  #Put results into dictionary
  output['left_inclusions']=left_inclusions
  output['right_inclusions']=right_inclusions
  left_ids="".join([',%d'%t[0] for t in left_inclusions])
  right_ids="".join([',%d'%t[0] for t in right_inclusions])
  output['left_ids']=left_ids
  output['right_ids']=right_ids
  output['all_inc_ids']=[t[0] for t in left_inclusions+right_inclusions]
  output['actual_number_inclusions']=len(output['all_inc_ids'])
  return output

def calc_tmpl_model(mesh_data, model_values):
  """Compute the values that go into the model template, from the relevant inputs"""
  output={}
  #Direct copy
  for k in ['all_inc_ids']:
    output[k]=mesh_data[k]
  for k in ['basename','meshname','modelname',
            'left_conc_A','left_conc_B',
            'right_conc_A','right_conc_B',
            'left_potential','right_potential',
            'r_inc_conc_A','r_inc_conc_B','r_inc_potential',
            'nr_inc_conc_A','nr_inc_conc_B','nr_inc_potential']:
    output[k]=model_values[k]
  #Calculations
  output['debye_length']=1.0/model_values['kappa']
  #Separate reactive and non-reactive inclusions
  #for now, assume only the first inclusion is reactive
  output['reactive_inc_ids']=[mesh_data['all_inc_ids'][0]]
  output['nonreactive_inc_ids']=mesh_data['all_inc_ids'][1:]
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

