"""For generating plots from simulation data"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp
import sys

#Site packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Local
import folderstructure as FS
import common
import simulator_general
import collect_results

class PlotSeries(common.ParameterSet):
  """Data for a single series on a plot

  Attributes:

    - xvals = array of x-values
    - yvals = array of y-values
    - label = legend label
    - metadata = other parameters needed to identify the data series
  """
  __slots__=['xvals','yvals','label','metadata']
  def add_to_axes(self,ax,fmt,**kwd):
    """Plot this series on the specified axes

    This is a wrapper for ax.plot

    Arguments:

      - ax = matplotlib Axes object
      - fmt = matplotlib format string
      - \**kwd = other keyword arguments for Axes.plot

    Returns:

      - The result of call to ax.plot"""
    return ax.plot(self.xvals,self.yvals,fmt,label=self.label,**kwd)

class PlotFigure(common.ParameterSet):
  """Data for a single matplotlib figure

  This is for a plot with a single set of axes.

  Attributes:

    To be read in from yaml file:

      - figsize = pair of numbers representing figure size, in inches: (Width, Height)
      - filename = name of the output file to be created, as string
      - prepfunctions = sequence of method calls used to generate additional data, etc.

        The available method names usually start with 'prep_'

      - plotfunctions = sequence of method calls used to generate plot

        The available method names usually start with 'plot_'

      - xlabel = x-axis label, as string
      - ylabel = y-axis label, as string
      - title = plot title, as string
      - fmts = list of format specifier strings

    To be created by methods:

      - datafiles = dictionary of loaded data files
      - outfpath = path to output file
      - series = sequence of PlotSeries instances
      - fig = matplotlib Figure for the generated plot
      - ax = matplotlib Axes for the plot
      - info = dictionary of miscellaneous data"""
  __slots__=['figsize','filename','prepfunctions','plotfunctions','xlabel','ylabel','title','fmts','outfpath','datafiles','series','fig','ax','info']
  _config_attrs=['figsize','filename','prepfunctions','plotfunctions','xlabel','ylabel','title','fmts']
  _outputfile_attrs=['outfpath']
  _taskname_src_attr='outfpath'
  
  def __init__(self,**kwd):
    #Initialization from base class
    super(PlotFigure, self).__init__(**kwd)
    #Find the input and output files
    self.locate_data()
    self.outfpath=osp.join(self.outdir(),self.filename)

  @property
  def _more_inputfiles(self):
    return list(self.datafiles.values())
  
  def execute_commandseq(self,attrname):
    """Execute the command sequence

    Arguments:

      - attrname = name of attribute containing the command sequence"""
    if getattr(self,attrname,None) is not None:
      for cmd in getattr(self,attrname,[]):
        #Function name and arguments
        funcname, kwargs = cmd
        #Call it
        try:
          getattr(self,funcname)(**kwargs)
        except Exception as einst:
          print("Excption occured for command: %s"%str(cmd), file=sys.stderr)
          raise einst
      
  def run(self):
    """Create the plot."""
    print(self.outfpath)
    
    #Load the data we need to generate the plot
    self.load_data()
    
    #Initialize the figure at the size requested
    self.fig = plt.figure(figsize=self.figsize)
    
    #Get the axes
    self.ax=self.fig.gca()
    
    #Call the preparation functions
    self.execute_commandseq('prepfunctions')
    
    #Add the available series to the axes
    self.plot_basic_series()
    
    #Call the requested plot functions
    self.execute_commandseq('plotfunctions')
    
    #Save the figure
    if not osp.isdir(self.outdir()):
      os.makedirs(self.outdir())
    self.fig.savefig(self.outfpath)
    
    #Close the figure
    plt.close(self.fig)
    
    #Done
    return

  def plot_basic_series(self):
    """A simple plot."""
    for i,sr in enumerate(self.series):
        o=sr.add_to_axes(self.ax,self.fmts[i])
    if getattr(self,'title',None) is not None:
      o=self.ax.set_title(self.title)
    if getattr(self,'xlabel',None) is not None:
      o=self.ax.set_xlabel(self.xlabel)
    if getattr(self,'ylabel',None) is not None:
      o=self.ax.set_ylabel(self.ylabel)
    return   

  def plot_axmethod(self,method,kwargs=None):
    """Call a method of the axes.

    Arguments:

      - method = name of Axes method to call, as string
      - kwargs = arguments dictionary for the method"""
    f=getattr(self.ax,method)
    if kwargs is None:
      kwargs = {}
    f(**kwargs)
    return

  def plot_hline(self,locspec,kwargs=None):
    """Add a horizontal line with a value from info

    Arguments:

      - locspec = sequence of keys in the info dictionary to locate the y-value
      - kwargs = keyword arguments for ax.axhline"""
    yval=common.nested_location(self.info,locspec)
    if kwargs is None:
      kwargs = {}
    self.ax.axhline(yval,**kwargs)
    return

  def plot_vline(self,locspec,kwargs=None):
    """Add a vertical line with a value from info

    Arguments:

      - locspec = sequence of keys in the info dictionary to locate the x-value
      - kwargs = keyword arguments for ax.axvline"""
    xval=common.nested_location(self.info,locspec)
    if kwargs is None:
      kwargs = {}
    self.ax.axvline(xval,**kwargs)
    return


class ModelPlotFigure(PlotFigure):
  """Data for a single model plot

  Attributes:

    To be read in from yaml file:

      - plotname = plot name in outdata file holding data series
      - modelname = name of model

    To be created by methods:

      (none)"""
  __slots__=['plotname','modelname']
  _config_attrs=PlotFigure._config_attrs+['plotname','modelname']

  def outdir(self):
    return osp.join(FS.postprocfolder,self.basename,self.modelname)

  def datadir(self):
    return osp.join(FS.solnfolder,self.basename,self.modelname)

  def locate_data(self):
    datadir=self.datadir()
    self.datafiles={'pklfile':osp.join(datadir,'outdata.pkl'), 'infofile':osp.join(datadir,FS.infofile)}

  def load_data(self):
    """Load the data for the plot."""
    
    #Load the data series
    outdata=simulator_general.OutData.from_pickle(self.datafiles['pklfile'])
    self.series=outdata.plots[self.plotname]
    
    #Load the info
    self.info=common.readyaml(self.datafiles['infofile'])
    
    return

class CollectionPlotFigure(PlotFigure):
  """Data for a single collection plot

  Attributes:

    To be read in from yaml file:

      - calcfunctions = sequence of calculation functions to be called before generating plot
      - seriesdefs = sequence of series definitions (xcol, ycol, label),
        where the columns specify the DataFrame columns containing values for the series.
        The label is optional.

    To be created by methods:

      - df = the DataFrame"""
  __slots__=['calcfunctions','seriesdefs','df']
  _config_attrs=PlotFigure._config_attrs+['calcfunctions','seriesdfs']

  def outdir(self):
    return osp.join(FS.postprocfolder,self.basename)

  def locate_data(self):
    self.datafiles={'dataframe': osp.join(FS.postprocfolder,self.basename,collect_results.collected_df_fname)}

  def load_data(self):
    """Load the data for the plot."""
    #Load the DataFrame
    self.df=pd.read_pickle(self.datafiles['dataframe'])
    
    #Initialize empty info
    self.info={}
    
    #Do the requested calculations to add new columns
    self.execute_commandseq('calcfunctions')
    
    #Add the requested columns in as series
    self.series=[]
    for sdef in self.seriesdefs:
      assert len(sdef) >= 2, "Inadequate series definition: %s"%str(sdef)
      if len(sdef)>3 and len(sdef[3])>0:
        qdf=self.df.query(sdef[3])
      else:
        qdf=self.df
      sdef_dict={'xvals':qdf[sdef[0]],'yvals':qdf[sdef[1]]}
      if len(sdef)>2:
        sdef_dict['label']=sdef[2]
      else:
        sdef_dict['label']=''
      self.series.append(PlotSeries(**sdef_dict))
    
    return

  def calc_Dratio(self):
    def calc_ratio(row):
      return row['Deff']/row['D_bulk']
    self.df['ratio_D']=self.df.apply(calc_ratio,axis=1)
    return

  def prep_series_equality(self):
    pdser=self.df['free_volume_frac']
    vals=[pdser.min(),pdser.max()]
    ser=PlotSeries(xvals=vals,yvals=vals,label="1:1")
    self.series.append(ser)
    return