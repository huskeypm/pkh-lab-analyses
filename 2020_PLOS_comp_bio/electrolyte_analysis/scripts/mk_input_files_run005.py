"""Create the input files from templates for run004

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
  - inc_conc_A, inc_conc_B: concentration values for Dirichlet boundary conditions on the inclusions
  - inc_potential: electric potential for Dirichlet boundary condition on the inclusions

Classes controlling fixed and random circular inclusions:
  - InclusionData: base class for the two sets of inclusions
  - LeftInclusions, RightInclusions: define parameters for the inclusions in each portion of the mesh
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

basename='run005'

#The number of characters in this string is the number of runs done in parallel
runlist='ABCDE'

#Parameter values: each combination of entries is used
mesh_vars=odict([
  ('randseed', range(5)),
  ('S', [12.0]),
  ('d', [8.0]),
  ('H', [40.0]),
  ('slot_length', [90.0]),
  ('mscale', [1.0]),
  ('mcar1', [5.0]),
  ('mcar2', [2.0]),
  ('mcar3', [1.0]),
  ('mcar4', [0.2]),
  ('incNum', range(20, 270, 20))
])

model_vars=odict([
  ('left_conc_A', [6.0e-4]),
  ('left_conc_B', [0.0]),
  ('right_conc_A', [6.0e-4]),
  ('right_conc_B', [0.0]),
  ('left_potential', [0.0]),
  ('right_potential', [0.0]),
  ('inc_conc_A', [0.0]),
  ('inc_conc_B', ['null']),
  ('inc_potential', [0.0, 25e-3, -25e-3]),
  ('kappa', [1.0])]
)

#  ('', []),

inc_id_start=1000
maxtries=1e4
class InclusionData(object):
  def __init__(self,mesh_values):
    mesh=Namespace(**mesh_values)
    self.mesh=mesh
  def gen_inclusions(self):
    global inc_id
    bbox=self.box
    centers,rvals=random_inclusions_2D.random_circles(self.fixed,self.nrandom,
                                                      (bbox.xmin,bbox.xmax),(bbox.ymin,bbox.ymax),self.radius_limits,
                                                      bbox,maxtries,self.eps)
    inclusions=[]
    for i,cen in enumerate(centers):
      inclusions.append((inc_id,cen[0],cen[1],rvals[i]))
      inc_id+=10
    return inclusions

def incCalc(incNum):
  lInc = incNum / 2
  rInc = lInc + (incNum % 2)
  return (lInc, rInc)

#incNum = 9
#lInc = incNum / 2
#rInc = lInc + (incNum % 2)

class LeftInclusions(InclusionData):
  def __init__(self,mesh_values):
    mesh=Namespace(**mesh_values)
    self.mesh=mesh
    self.box=random_inclusions_2D.Rect(-mesh.slot_length/2,-mesh.d/2,0,mesh.d/2)
    self.fixed=[(-mesh.slot_length/2 + 1.5,0,1)]
    self.nrandom=int(incCalc(mesh_values['incNum'])[0] - len(self.fixed))
    self.radius_limits=(0.5,0.5)
    self.eps=0.25

class RightInclusions(InclusionData):
  def __init__(self,mesh_values):
    mesh=Namespace(**mesh_values)
    self.mesh=mesh
    self.box=random_inclusions_2D.Rect(0,-mesh.d/2,mesh.slot_length/2,mesh.d/2)
    self.fixed=[]
    self.nrandom=int(incCalc(mesh_values['incNum'])[1])
    self.radius_limits=(0.5,0.5)
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

#Loop over meshes
modelid = 0
meshvals_iter=product(*list(mesh_vars.values()))
for meshid,meshvals in enumerate(meshvals_iter):
  meshname='%s_mesh_%03d'%(basename,meshid)
  mesh_values=odict(zip(mesh_vars.keys(),meshvals))
  #Generate inclusions
  random_inclusions_2D.reseed(mesh_values['randseed'])
  inc_id=inc_id_start
  incdata_l=LeftInclusions(mesh_values)
  incdata_r=RightInclusions(mesh_values)
  left_inclusions=incdata_l.gen_inclusions()
  right_inclusions=incdata_r.gen_inclusions()
  mesh_values['left_inclusions']=left_inclusions
  mesh_values['right_inclusions']=right_inclusions
  left_ids="".join([',%d'%t[0] for t in left_inclusions])
  right_ids="".join([',%d'%t[0] for t in right_inclusions])
  mesh_values['left_ids']=left_ids
  mesh_values['right_ids']=right_ids
  all_inc_ids=[t[0] for t in left_inclusions+right_inclusions]
  #Create the file from the template
  geo_outfpath=osp.join(geo_outdir,'%s.geo'%meshname)
  apply_template(geo_tmpl,geo_outfpath,mesh_values)
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
    model_values['debye_length']=1.0/model_values['kappa']
    model_values['all_inc_ids']=all_inc_ids
    #Add model to yaml document
    yaml_doc=model_tmpl.render(**model_values)
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

