"""Defining the data series that go on plots"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from copy import deepcopy

#Site packages
import numpy as np

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.commandseq import CommandSequenceRequest

class PlotSeries(object):
  """Data for a single series on a plot

  Attributes:

    - xvals = array of x-values, required
    - yvals = array of y-values, required
    - label = legend label, optional
    - metadata = other parameters needed to identify the data series, optional
    - errors = error values in y, optional
  """
  _expected_attrs=['xvals','yvals','label','metadata','errors']
  def __init__(self,xvals,yvals,label=None,metadata=None,errors=None):
    for attrname in self._expected_attrs:
      setattr(self,attrname,locals()[attrname])
  def add_to_axes(self,ax,fmt,newlabel=None,**kwd):
    """Plot this series on the specified axes

    This is a wrapper for ax.plot

    Arguments:

      - ax = matplotlib Axes object
      - fmt = matplotlib format string
      - newlabel = label to use for the legend instead of the series-provided one
      - \**kwd = other keyword arguments for Axes.plot

    Returns:

      - The result of call to ax.plot"""
    #Handle label
    if newlabel is None:
      label=getattr(self,'label',None)
      if label is None:
        label=''
    else:
      label=newlabel
    #Use the right plot method
    if getattr(self,'errors',None) is not None:
      res = ax.errorbar(self.xvals,self.yvals,self.errors,fmt=fmt,label=label,**kwd)
    else:
      res = ax.plot(self.xvals,self.yvals,fmt,label=label,**kwd) 
    return res
  def __getstate__(self):
    """Used for pickling, and converting to yaml"""
    state={}
    for attrname in self._expected_attrs:
      state[attrname]=getattr(self,attrname)
    return state
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)

#Register for loading from yaml
yaml_manager.register_classes([PlotSeries])
