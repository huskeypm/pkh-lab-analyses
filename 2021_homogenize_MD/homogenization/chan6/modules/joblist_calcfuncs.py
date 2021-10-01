"""Calculation functions for the job list"""

def calc_resfile(self,row,setname,fname):
  resfile="solutions/{setname}/{tag}/{job_id}/{fname}".format(setname=setname,fname=fname,**row)
  return resfile

#List of functions to be bound as methods
request_methods=[calc_resfile]