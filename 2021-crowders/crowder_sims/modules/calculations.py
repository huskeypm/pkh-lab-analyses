"""Post-processing calculations"""

import numpy as np

from simproc.postproc.plotseries import PlotSeries

def calc_ratio_column(self, numcol, dencol, newcol, dfpath="df", const=1.0):
  """Compute the porosity column of the dataframe

  Arguments:
  
    - numcol = column name for the numerator
    - dencol = column name for the denominator
    - newcol = name of new column storing results, as string
    - dfpath = optional, attribute path to the dataframe
    - const = optional, constant value to multiply the ratio by

  New column added to the dataframe.
  No return value.
  No output files."""
  df=self.get_nested(dfpath)
  def divide(row):
    return const*row[numcol]/row[dencol]
  df[newcol]=df.apply(divide,axis=1)
  return

#List of functions to be bound as methods
request_methods=[calc_ratio_column]
