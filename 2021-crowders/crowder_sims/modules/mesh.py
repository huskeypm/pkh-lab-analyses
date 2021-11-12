"Customization of mesh"

#Site packages
import pandas as pd

#From simproc
from simproc.requesthandler.nested import WithNested

#This directory
import lattice_inclusions

#Constants
mesh_direct_fields=['pore_radius','pore_length','RevW','RevH','RevD',
                    'mscale','mcar1','mcar2','mcar3','mcar4','inclusionlist']

inc_id_step = 100

inclusion_columns=['ID','centerX','centerY','centerZ','radius']

def mesh_get_template_input(self):
  """Compute the values that go into the mesh template, from the relevant inputs
  
  self.data is expected to contain the following:

    - all fields in ``mesh_direct_fields``
    - optionally, a dictionary ``lattice`` to define the lattice inclusions

  """
  output={}
  #Split static inclusions into parts above and below pore center
  lower_inclusions=[t for t in self.data['inclusionlist'] if t[3] < 0.0]
  upper_inclusions=[t for t in self.data['inclusionlist'] if t[3] > 0.0]
  #Compute lattice inclusions
  if 'lattice' in self.data.keys():
    #Get necessary parameters
    data=WithNested(**self.data)
    data.lattice=WithNested(**data.lattice)
    incrad=data.lattice.incrad
    idstart=data.lattice.id_start
    cylrad=data.pore_radius
    buffer=data.lattice.sep_buffer
    lowend=(0.0,0.0,-data.pore_length/2.0)
    highend=(0.0,0.0,data.pore_length/2.0)
    origin=(0.0,0.0,0.0)
    perturbation_rounds=data.lattice.get('perturbation_rounds',0)
    perturbation_stdev=data.lattice.get('perturbation_stdev',0)
    #Generate inclusions in lower half of pore
    new_lower=lattice_inclusions.generate_lattice_inclusions_in_cylinder(
        idstart,incrad,cylrad,buffer,lowend,origin,lower_inclusions,inc_id_step,
        perturbation_rounds,perturbation_stdev)
    #Generate inclusions in upper half of pore
    lastid=new_lower[-1][0]
    idstart=lastid+inc_id_step
    new_upper=lattice_inclusions.generate_lattice_inclusions_in_cylinder(
        idstart,incrad,cylrad,buffer,origin,highend,upper_inclusions,inc_id_step,
        perturbation_rounds,perturbation_stdev)
    #Update all inclusion lists
    self.data['inclusionlist']+=new_lower+new_upper
    upper_inclusions+=new_upper
    lower_inclusions+=new_lower
    self.inclusions_df=pd.DataFrame(self.data['inclusionlist'],columns=inclusion_columns)
  #Direct copy
  for k in mesh_direct_fields:
    output[k]=self.data[k]
  #Inclusion ID lists (some are in string form)
  output['top_inclusion_ids']="".join([',%d'%t[0] for t in upper_inclusions])
  output['bot_inclusion_ids']="".join([',%d'%t[0] for t in lower_inclusions])
  # output['all_inc_ids']=[t[0] for t in self.data['inclusionlist']]
  return output

#List of functions to be bound as methods
request_methods=[mesh_get_template_input]
