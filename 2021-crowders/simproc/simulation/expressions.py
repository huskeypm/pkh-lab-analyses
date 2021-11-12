"""Support module for fenics: UserExpression subclasses"""

#Site packages
import fenics as fem


class VaryingScalarByCell(fem.UserExpression):
  def __init__(self,meshfunc,prop_map,**kwargs):
    self.meshfunc=meshfunc
    self.prop_map=prop_map
    super().__init__(**kwargs)
  def value_shape(self):
    return ()
  def eval_cell(self, values, x, cell):
    key=self.meshfunc[cell.index]
    values[0]=self.prop_map[key]

class VaryingVectorByCell(fem.UserExpression):
  def __init__(self,meshfunc,prop_map,spatial_dims,**kwargs):
    self.meshfunc=meshfunc
    self.prop_map=prop_map
    self.spatial_dims=spatial_dims
    super().__init__(**kwargs)
  def value_shape(self):
    return (self.spatial_dims,)
  def eval_cell(self, values, x, cell):
    key=self.meshfunc[cell.index]
    for i in range(self.spatial_dims):
        values[i]=self.prop_map[key][i]

class VaryingMatrixByCell(fem.UserExpression):
  def __init__(self,meshfunc,prop_map,spatial_dims,**kwargs):
    self.meshfunc=meshfunc
    self.prop_map=prop_map
    self.spatial_dims=spatial_dims
    super().__init__(**kwargs)
  def value_shape(self):
    return (self.spatial_dims,self.spatial_dims)
  def eval_cell(self, values, x, cell):
    key=self.meshfunc[cell.index]
    for i in range(self.spatial_dims):
      for j in range(self.spatial_dims):
        values[self.spatial_dims*j+i]=self.prop_map[key][i][j]
