"""Create the input files from templates for run001

Variables that must be set in the mesh template file (.geo):

  - inclusion_dist: distance from pore surface to inclusion

Variables that must be set in the model template file (.yaml):

  - meshname
  - modelname
  - sphere_conc: species concentration at the sphere surface
  - top_conc: species concentration on the "top" reservoir surface
  - bot_conc: species concentration on the "bottom" reservoir surface
  - central_pore_pot: electric potential on the surface of the central pore
  - sphere_pot: electric potential on the surface of the sphere
  - debye_length: debye length for the electrostatic potential (1/kappa)
  """

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from copy import deepcopy
import os
import os.path as osp
import sys

#Site packages
import yaml

#Locate source folders
if 'SRCLOC' in os.environ.keys():
  srcfolder=osp.normpath(osp.abspath(os.environ['SRCLOC']))
else:
  srcfolder=osp.abspath(osp.split(__file__)[0])

#add python code folder(s) to path
if not srcfolder in sys.path:
  sys.path.append(srcfolder)

import folderstructure as FS

basename='run001'

#Constant model parameters
model_consts={'sphere_conc': 0.0,
              'top_conc': 6.0e-4,
              'bot_conc': 0.0,
              'sphere_pot': 0.0}

geo_template_file = osp.join(FS.datafolder,"mesh/templates/%s.geo.fmt"%basename)
mesh_yaml_outfile=osp.join(FS.params_mesh_folder,'%s.yaml'%basename)
geo_outdir=osp.join(FS.geofolder,basename)
meta_outdir=osp.join(FS.meshmeta_outfolder,basename)
for dirname in [geo_outdir, meta_outdir]:
  if not osp.isdir(dirname):
    os.mkdir(dirname)
model_template_file = osp.join(FS.datafolder,"model/templates/model_%s.yaml.fmt"%basename)
model_yaml_outfile=osp.join(FS.params_model_folder,'%s.yaml'%basename)
codebook_outfile='codebook_%s.yaml'%basename

dist_list = [0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
pore_pot_list = [0.0, 25.0e-3,-25.0e-3]
kappa_list=[0.01, 0.1, 1.0, 10.0, 100.0]

def apply_template(tmpl_data,out_fpath,values):
  """Apply the given data to the specified template to produce an output file.
  
  Arguments:
  
    - tmpl_data = template data, as string set up for use with the .format method
    - out_fpath = path to output file, as string
    - values = dictionary of values to put into the template, {template variable name: value}"""
  #Substitute output
  out_data=tmpl_data.format(**values)
  #Write output file
  with open(out_fpath,'w') as fh:
    fh.write(out_data)

#Read the yaml input templates
with open(geo_template_file,'r') as fh:
  geo_tmpl_data=fh.read()
with open(model_template_file,'r') as fh:
  model_tmpl_data=fh.read()

#Initialize yaml output
mesh_yaml_out="""%YAML 1.2\n"""
model_yaml_out="""%YAML 1.2\n"""
codebook_out="""%YAML 1.2\n---\n"""

modelid=0
#Loop over meshes
#Only one varying parameter: inclusion distance
for meshid,dist in enumerate(dist_list):
  #Create the file from the template
  mesh_values={'inclusion_dist':dist}
  meshname='%s_mesh_%03d'%(basename,meshid)
  geo_outfpath=osp.join(geo_outdir,'%s.geo'%meshname)
  apply_template(geo_tmpl_data,geo_outfpath,mesh_values)
  #Create a basic mesh metadata file
  meta_str='\n'.join(['%s: %s'%tup for tup in mesh_values.items()])
  with open(osp.join(meta_outdir,'%s.yaml'%meshname),'w') as fh:
    fh.write(meta_str)
  #Add to the yaml output
  yaml_doc="---\nmeshname: %s\ngeomdefname: null\ntmplvalues:\n"%(meshname)
  yaml_doc+="\n".join(["  %s: %s"%tup for tup in mesh_values.items()])+"\n"
  mesh_yaml_out+=yaml_doc
  #Loop over simulations to be done with this mesh
  for pore_pot in pore_pot_list: #Each pore potential
    #A potential of zero means no dependence on kappa
    if pore_pot == 0.0:
      kappa_iter=[1]
    else:
      kappa_iter=kappa_list
    for kappa in kappa_iter: #Each kappa value
      modelname='%s_model_%03d'%(basename,modelid)
      #Add this to the codebook
      codebook_out+="%s:\n  inc_dist: %s\n  pore_pot: %s\n  kappa: %s\n"%(modelname,dist,pore_pot,kappa)
      #Get the dictionary of values for the template
      model_values=deepcopy(model_consts)
      model_values['meshname']=meshname
      model_values['modelname']=modelname
      model_values['central_pore_pot']=pore_pot
      model_values['debye_length']=1.0/kappa
      #Add model to yaml document
      yaml_doc=model_tmpl_data.format(**model_values)
      model_yaml_out+=yaml_doc
      #Next model
      modelid+=1


#Write the output yaml files
with open(codebook_outfile,"w") as fh:
  fh.write(codebook_out)
with open(mesh_yaml_outfile,"w") as fh:
  fh.write(mesh_yaml_out)
with open(model_yaml_outfile,'w') as fh:
  fh.write(model_yaml_out)


