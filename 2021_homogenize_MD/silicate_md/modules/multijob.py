"Preparation of multiple jobs with varying methane location"

#Standard library
import itertools

#Site packages
import numpy as np
import pandas as pd

#From simproc

#This directory

#Constants
direct_fields_required=[
  "setname","stage","cluster","restrain","nsteps_min",
  "nsteps_equil","nsteps_prod","timestep","spring","temp","boxspace",
  "has_silicate","has_methane","silicate_id",
  "do_basic_production","extension_A","extension_B"]
direct_fields_optional=["xout_steps","xout_grps","postproc_template","postproc_template_list",
  "max_duration", "node_exclusions"]
joblist_columns=['jobname','atm_x','atm_y','atm_z']
coordkeys="xyz"
linedata_keys=["hist_min_y","hist_max_y","num_bins"]
planedata_keys=["hist_min_x","hist_max_x","num_bins_x","hist_min_y","hist_max_y","num_bins_y"]

def get_template_input(self):
  """Compute the values that go into the template, from the relevant inputs
  
  self.data is expected to contain the following:

    - all fields in ``direct_fields_required``
    - optionally, any fields in ``direct_fields_optional``
    - start: dictionary of starting coordinate values, with keys x, y, and z
    - stop: dictionary of stopping coordinate values, with keys x, y, and z
        (depending on the step size, the stopping coordinate may be exceeded,
        or may not be reached.)
    - step: dictionary of step size values, with keys x, y, and z
        (for any direction with start=stop, the step may be omitted)
  """
  output={}
  #Direct copy
  for k in direct_fields_required:
    assert k in self.data, "Template data for %s missing required field %s"%(getattr(self,"name","(unnamed)"),k)
    output[k]=self.data[k]
  for k in direct_fields_optional:
    if k in self.data:
      output[k]=self.data[k]
  #Additional validation
  if not "postproc_template_list" in self.data.keys():
    assert "postproc_template" in self.data.keys(), "Must define either postproc_template or postproc_template_list"
  #Check the validity of the start, stop, and step values
  for k in coordkeys:
    assert k in self.data["start"], "'start' is missing %s"%k
    assert k in self.data["stop"], "'stop' is missing %s"%k
    assert self.data["stop"][k]>=self.data["start"][k], "'stop' is less than 'start' for %s"%k
    if self.data["stop"][k]>self.data["start"][k]:
      assert k in self.data["step"], "'step' is missing %s"%k
  #Generate all the coordinate values
  self.coordvals={}
  for k in coordkeys:
    startval=self.data["start"][k]
    stopval=self.data["stop"][k]
    if startval==stopval:
      vlist=[startval]
    else:
      stepval=self.data["step"][k]
      nsteps=int(np.ceil((stopval-startval)/stepval))
      vlist=[startval + idx * stepval for idx in range(nsteps)]
    self.coordvals[k]=vlist
  #Generate the list of jobs
  joblist=[]
  jobnum=0
  iterlist=[self.coordvals[k] for k in coordkeys]
  for xval,yval,zval in itertools.product(*iterlist):
    jobnum += 1
    jobname = "%04d"%jobnum
    jobtup=[jobname,xval,yval,zval]
    joblist.append(jobtup)
  #Store results as needed
  output["joblist"]=joblist
  self.joblist=joblist
  self.joblist_df=pd.DataFrame(joblist,columns=joblist_columns)
  #Preparations for 1D wham
  if 'linedata' in self.data.keys():
    #Generate list of windows for 1D wham
    self.windows={}
    lineno=0
    for xval,yval in itertools.product(self.coordvals["x"],self.coordvals["y"]):
      lineno+=1
      lineid="%03d"%lineno
      subset=self.joblist_df.loc[self.joblist_df['atm_x']==xval]
      subset=subset.loc[subset['atm_y']==yval]
      self.windows[lineid]=[[row.jobname,row.atm_z] for row in subset.itertuples()]
    #Generate list of lines for 1D wham
    self.linelist=[]
    line_data=[self.data['linedata'][k] for k in linedata_keys]
    for lineid in self.windows.keys():
      self.linelist.append([lineid]+line_data)
    output['linelist']=self.linelist
  #Preparations for 2D wham
  if 'planedata' in self.data.keys():
    #Generate list of windows for 2D wham
    self.windows2d={}
    planeno=0
    for xval in self.coordvals["x"]:
      planeno+=1
      planeid="%02d"%planeno
      subset=self.joblist_df.loc[self.joblist_df['atm_x']==xval]
      self.windows2d[planeid]=[[row.jobname,row.atm_y,row.atm_z] for row in subset.itertuples()]
    #Generate list of planes for 2D wham
    self.planelist=[]
    plane_data=[self.data['planedata'][k] for k in planedata_keys]
    for planeid in self.windows2d.keys():
      self.planelist.append([planeid]+plane_data)
    output['planelist']=self.planelist
  #Done
  return output

#List of functions to be bound as methods
request_methods=[get_template_input]
