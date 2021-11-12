"Customizations of model"

def assign_crowder_potential(self,potential,df_attrpath,idcol="ID",bc_attrpath="conditions_processed.dirichlet"):
  """Add the crowder inclusions to the Dirichlet boundary conditions list

  Arguments:

    - potential = value of the potential to use for the crowders, as a float, None to apply no Dirichlet condition
    - df_attrpath = attribute path to the dataframe listing the inclusions
    - idcol = name of dataframe column containing inclusion ID
    - bc_attrpath = attribute path for the Dirichlet boundary conditions dictionary
  """
  #Do not apply Dirichlet conditions if potential is None
  if potential is not None:
    #Get the inclusion list dataframe
    df=self.get_nested(df_attrpath)
    #Set up a dictionary of all the boundary conditions to add
    newdict={}
    for tup in df.itertuples(index=False,name='inclusions'):
      newdict[getattr(tup,idcol)]=potential
    #Append the new data to the Dirichlet boundary conditions dictionary
    #The inclusion list loaded above has the non-crowders in it as well,
    # so overwrite those entries with the potentials specified in the request,
    # rather than the crowder potential.
    bcdict=self.get_nested(bc_attrpath)
    newdict.update(bcdict)
    self.set_nested(bc_attrpath,newdict)
  #Done
  return

#List of functions to be bound as methods
request_methods=[assign_crowder_potential]
