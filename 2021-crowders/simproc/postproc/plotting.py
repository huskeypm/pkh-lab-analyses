"""Generate plots"""

#Standard library
import os

#Site packages
import numpy as np
import matplotlib as mpl
#Set matplotlib backend that won't require Tk
if mpl.get_backend() == 'TkAgg':
  mpl.use("Agg") #this line must appear before the following, and you can't change it later
import matplotlib.pyplot as plt
import pandas as pd

#This package
from ..requesthandler.commandseq import WithCommandsRequest
from ..requesthandler import schema
from ..requesthandler import yaml_manager
from .plotseries import PlotSeries
from ..requesthandler import logging

logger=logging.getLogger(__name__)

_FigureProperties_props_schema_yaml="""#FigureProperties
figattr: {type: attrpath}
figsize:
  type: array
  minItems: 2
  maxItems: 2
  items: {type: number}
dpi: {type: number}
facecolor: {}
edgecolor: {}
linewidth: {type: number}
frameon: {type: boolean}
subplotpars: {}
tight_layout: {type: boolean}
constrained_layout: {type: boolean}
outfpath: {type: pathlike}
savekwargs: {type: object}
"""

class FigureProperties(schema.SelfValidating):
  """Properties of a figure
  
  Arguments to ``plt.figure`` are allowed as attributes.
  
  Other attributes:
  
    - figattr = attribute path for storing result, defaults to "fig"
    - outfpath = file path to save to once complete
    - savekwargs = keyword arguments dictionary to pass to ``fig.savefig``"""
  _validation_schema=schema.SelfValidating.update_schema(_FigureProperties_props_schema_yaml)
  _validation_schema.required=['outfpath']
  _matplotlib_props=['figsize','dpi','facecolor','edgecolor','linewidth','frameon','subplotpars','tight_layout','constrained_layout']

schema.extra_types_dict['FigureProperties']=(FigureProperties,)

_AxesProperties_props_schema_yaml="""#AxesProperties
axattr: {type: attrpath}
figattr: {type: attrpath}
figlist:
  type: array
  items: {type: attrpath}
nrows: {type: integer}
ncols: {type: integer}
index: {type: integer}
title:
  anyOf:
    - {type: string}
    - {type: object}
xlabel:
  anyOf:
    - {type: string}
    - {type: object}
ylabel:
  anyOf:
    - {type: string}
    - {type: object}
xlim:
  type: array
  minItems: 2
  maxItems: 2
  items: {type: number}
ylim:
  type: array
  minItems: 2
  maxItems: 2
  items: {type: number}
kwargs: {type: object}
"""

class AxesProperties(schema.SelfValidating):
  """Properties of a set of axes (a single plot in a potentially multiplot figure)
  
  Properties:
  
    - axattr = attribute path for storing result, defaults to "ax"
    - figattr = attribute path for figure to put the axes on (do not combine with ``figlist`` property)
    - figlist = list of attribute paths for figures to put the axes on (do not combine with ``figattr`` property)
    - nrows, ncols, index = integers specifying subplot location (all default to 1)
    - title = plot title, as string
    - xlabel = label for x-axis, as string OR dictionary of keyword arguments to ``ax.set_xlabel``
    - ylabel = label for y-axis, as string OR dictionary of keyword arguments to ``ax.set_ylabel``
    - xlim = pair of numbers to set as x-axis limits
    - ylim = pair of numbers to set as y-axis limits
    - kwargs = keyword arguments to pass to plt.figure.subplots"""
  _validation_schema=schema.SelfValidating.update_schema(_AxesProperties_props_schema_yaml)
  _validation_schema.required=[]

schema.extra_types_dict['AxesProperties']=(AxesProperties,)

_SeriesProperties_props_schema_yaml="""#SeriesProperties
seriesattr: {type: attrpath}
axattr: {type: attrpath}
axlist:
  type: array
  items: {type: attrpath}
fmt: {type: string}
newlabel: {type: string}
kwargs: {type: object}
"""

class SeriesProperties(schema.SelfValidating):
  """Properties of a data series to put on a set of axes
  
  Properties:
  
    - seriesattr = nested path to the series instance to add to the axes
    - axattr = attribute path for axes to put the series on (do not combine with ``axlist`` property)
    - axlist = list of attribute paths for figures to put the axes on (do not combine with ``axattr`` property)
    - fmt = matplotlib format string
    - newlabel = label to use for the legend instead of the series-provided one
    - kwargs = other keyword arguments for Axes.plot"""
  _direct_attrs=['fmt','newlabel']
  _validation_schema=schema.SelfValidating.update_schema(_SeriesProperties_props_schema_yaml)
  _validation_schema.required=['seriesattr']

schema.extra_types_dict['SeriesProperties']=(SeriesProperties,)

_FigureRequest_props_schema_yaml="""#FigureRequest
loadfiles:
  type: object
  additionalProperties:
    items: pathlike
prepcommands: {type: array}
rcparams: {type: object}
figures:
  type: array
  items: {type: FigureProperties}
axes:
  type: array
  items: {type: AxesProperties}
series:
  type: array
  items: {type: SeriesProperties}
basecommands: {type: array}
plotcommands: {type: array}
_all_axattr: {type: array}"""


class FigureRequest(WithCommandsRequest):
  """For generating matplotlib figures
  
  User-defined attributes:
  
    - loadfiles = dictionary of pickle files to load (yaml files may be loaded using a prepcommand), usually to load series for plotting
        The dictionary is formatted as {attribute: file path or locator},
        where the attribute to store the loaded data could also be an attribute path.
    - prepcommands = sequence of commands to execute prior to plotting, such as for defining series (see below)
    - rcparams = dictionary of matplotlib rcParams to be set
    - figures = dictionary of figure definitions {attribute path: FigureProperties instance}
    - axes = dictionary of axes definitions, {attribute path: AxesProperties instance}
    - series = list of series definitions, each a SeriesProperties instance
    - basecommands = sequence of commands to execute before series have been added to the axes (useful for contour plots)
    - plotcommands = sequence of commands to execute after series have been added to the axes

    prepcommands, basecommands, and plotcommands are all command lists.
    Each command in the list is a pair (cmdname, arguments), where:

      - cmdname = name of the object's method to call, as a string
      - arguments = dictionary of arguments to the method: {argname: value,...}
  
  Calculated attributes:

    - _all_axattr = list of all attributes containing matplotlib Axes objects added to the figure
  """
  _self_task=True
  _config_attrs=('loadfiles','prepcommands','basecommands','plotcommands','rcparams','figures','axes','series')
  _validation_schema=WithCommandsRequest.update_schema(_FigureRequest_props_schema_yaml)
  _validation_schema.required=[]
  def __init__(self,**kwargs):
    #Initialization from base class
    super(FigureRequest, self).__init__(**kwargs)
    #Get input and output files from the command sequences
    self.init_command_sequence('prepcommands')
    self.init_command_sequence('plotcommands')
    #Get input files from the loadfiles attribute
    self._more_inputfiles+=[fp for k,fp in getattr(self,'loadfiles',{}).items()]
    #Get output files from the figures attribute
    self._more_outputfiles+=[figprops.outfpath for figprops in getattr(self,'figures',[])]
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Load data files
    for attrpath,floc in getattr(self,'loadfiles',{}).items():
      self.load_pickle(floc,attrpath)
    #Run prepcommands
    self.process_command_sequence(attrpath='prepcommands',singlefunc=None,positional=False)
    #Set rcParams
    self.set_rcparams(**getattr(self,'rcparams',{}))
    #Set up figures
    for figprops in getattr(self,'figures',[]):
      self.figure_figprops(figprops)
    #Set up axes
    self._all_axattr=[]
    for axprops in getattr(self,'axes',[]):
      self.axes_axprops(axprops)
    #Run basecommands
    self.process_command_sequence(attrpath='basecommands',singlefunc=None,positional=False)
    #Add series to axes
    for serprops in getattr(self,'series',[]):
      self.add_seriesprops(serprops)
    #Run plotcommands
    self.process_command_sequence(attrpath='plotcommands',singlefunc=None,positional=False)
    #Save and close figures
    for figprops in getattr(self,'figures',[]):
      self.saveclose_figprops(figprops)
    #Done
    return
  def series_equalityline(self,outattr,span,label=None,metadata=None):
    """Generate a series for a 1:1 line
    
    Arguments:
    
      - outattr = nested path for storing the resulting series
      - span = pair of numbers for the lower and upper ends of the line
      - label = optional series label
      - metadata = optional series metadata"""
    series=PlotSeries(xvals=span,yvals=span,label=label,metadata=metadata)
    self.set_nested(outattr,series)
  def series_from_dataframe(self,dfpath,xcol,ycol,outattr,query=None,label=None,metadata=None,errcol=None):
    """Generate a series based on previously-loaded pandas dataframe

    Arguments:
    
      - dfpath = nested path to the DataFrame (NOT a file path or locator)
      - xcol = name of the DataFrame column to use for the x-values of the series
      - ycol = name of the DataFrame column to use for the y-values of the series
      - outattr = nested path for storing the resulting series
      - query = optional, DataFrame query to use to select only certain rows of the dataframe, as an argument to pd.DataFrame.query
      - label = optional series label
      - metadata = optional series metadata
      - errcol = optional name of DataFrame column to use for the y-error values of the series (symmetric)
    
    If you want to save the result to a pickle file, use a subsequent call to ``save_pickle``."""
    #Get the dataframe
    df = self.get_nested(dfpath)
    #Select particular rows if requested
    if query is None or len(query)==0:
      qdf=df
    else:
      qdf=df.query(query)
    #Get the data for the series from the dataframe
    xvals=qdf[xcol]
    yvals=qdf[ycol]
    if errcol is not None:
      errvals=qdf[errcol]
    else:
      errvals = None
    #Instantiate the series
    series=PlotSeries(xvals=xvals,yvals=yvals,label=label,metadata=metadata,errors=errvals)
    #Store result
    self.set_nested(outattr,series)
  def series_from_normalization(self,inpath,outattr,normpath,label=None,metadata=None):
    """Generate a PlotSeries by normalizing another one

    Arguments:

      - inattr = nested path to the input data series
      - outattr = nested path to the output data series
      - normpath = nested path to the normalization constant
      - label = optional series label
      - metadata = optional series metadata (added to the metadata from the other series)"""
    #Data from other series
    other=self.get_nested(inattr)
    xvals=other.xvals
    orig_y=np.array(other.yvals)
    out_meta=deepcopy(other.metadata)
    #Do the normalization
    factor=self.get_nested(normpath)
    yvals=orig_y/factor
    #Update metadata
    out_meta.update(metadata)
    #Instantiate the series
    series=PlotSeries(xvals=vals,yvals=yvals,label=label,metadata=out_meta)
    #Store result
    self.set_nested(outattr,series)
  def series_from_constant(self,inattr,outattr,span,label=None,metadata=None):
    """Generate a PlotSeries for a constant value
    
    Arguments:
    
      - inattr = nested path to the constant
      - outattr = nested path to the output data series
      - span = pair of numbers for the lower and upper ends of the line
      - label = optional series label
      - metadata = optional series metadata"""
    constval=self.get_nested(inattr)
    yvals=[constval, constval]
    series=PlotSeries(xvals=span,yvals=yvals,label=label,metadata=metadata)
    self.set_nested(outattr,series)
    return
  def series_closedbox(self,xvals,yvals,outattr,label=None,metadata=None):
    """Generate a PlotSeries making a closed box aligned with the axes

    Arguments:

      - xvals = min and max x values (or ``Stored`` pointing to them)
      - yvals = min and max y values (or ``Stored`` pointing to them)
      - outattr = nested path to the output data series
      - label = optional series label
      - metadata = optional series metadata"""
    xvals=self.get_stored(xvals)
    yvals=self.get_stored(yvals)
    xdata=[xvals[idx] for idx in [0,1,1,0,0]]
    ydata=[yvals[idx] for idx in [0,0,1,1,0]]
    series=PlotSeries(xvals=xdata,yvals=ydata,label=label,metadata=metadata)
    self.set_nested(outattr,series)
  def set_rcparams(self,**kwargs):
    """To set matplotlib rcParams"""
    for k,v in kwargs.items():
      mpl.rcParams[k]=v
    return
  def figure(self,figattr="fig",**kwargs):
    """To create a new matplotlib.figure instance

    Arguments:

      - figattr = nested path to store generated figure, defaults to "fig"
      - \**kwargs = keyword arguments to pass to plt.figure"""
    fig=plt.figure(**kwargs)
    self.set_nested(figattr,fig)
    return
  def figure_figprops(self,figprops):
    """Create a new figure from a FigureProperties instance"""
    kwargs={}
    for attr in figprops._matplotlib_props:
      if hasattr(figprops,attr):
        kwargs[attr]=getattr(figprops,attr)
    figattr=getattr(figprops,'figattr',"fig")
    self.figure(figattr,**kwargs)
    return
  def savefig(self,outfpath,figattr="fig",**kwargs):
    """Save the figure to file"""
    fig=self.get_nested(figattr)
    fig.savefig(self.renderstr(outfpath),**kwargs)
  def closefig(self,figattr="fig"):
    """Close the requested figure"""
    fig=self.get_nested(figattr)
    plt.close(fig)
  def saveclose_figprops(self,figprops):
    """Save and close the figure identified by its FigureProperties"""
    figattr=getattr(figprops,'figattr',"fig")
    fig=self.get_nested(figattr)
    savekwargs=getattr(figprops,'savekwargs',{})
    fig.savefig(self.renderstr(figprops.outfpath),**savekwargs)
    plt.close(fig)
  def figmethod(self,method,figattr="fig",outattr=None,**kwargs):
    """Call the given matplotlib Figure method

    Arguments:

      - method = method name
      - figattr = nested path to store generated figure, defaults to "fig"
      - \**kwargs = keyword arguments to the method"""
    fig=self.get_nested(figattr)
    f=getattr(fig,method)
    out=f(**kwargs)
    if outattr is not None:
      self.set_nested(outattr,out)
    return
  def axes(self,nrows=1,ncols=1,index=1,figattr="fig",axattr="ax",**kwargs):
    """To create a new matplotlib Axes instance

    Arguments:

      - nrows = number of rows of subplots, defaults to 1
      - ncols = number of columns of subplots, defaults to 1
      - index = index number of subplot, identifying which row and column (see matplotlib docs for explanation), defaults to 1
      - figattr = nested path to figure instance, defaults to "fig"
      - axattr = nested path to resulting Axes instance, defaults to "ax"
      - \**kwargs = keyword arguments to pass to plt.figure.subplots"""
    ax=self.get_nested(figattr).add_subplot(nrows,ncols,index,**kwargs)
    self.set_nested(axattr,ax)
    self._all_axattr.append(axattr)
    return
  def axes_axprops(self,axprops):
    """Create new axes from an AxesProperties instance"""
    #Axes attribute path
    axattr=getattr(axprops,'axattr',"ax")
    self._all_axattr.append(axattr)
    #Figure attribute paths
    if hasattr(axprops,'figattr'):
      assert not hasattr(axprops,'figlist'), "AxesProperties instance has both figattr (%s) and figlist (%s)"%(str(axprops.figattr),str(axprops.figlist))
      figlist=[axprops.figattr]
    else:
      if hasattr(axprops,'figlist'):
        figlist=axprops.figlist
      else:
        figlist=["fig"]
    #Subplot position
    nrows=getattr(axprops,'nrows',1)
    ncols=getattr(axprops,'ncols',1)
    index=getattr(axprops,'index',1)
    #Keyword arguments
    kwargs=getattr(axprops,'kwargs',{})
    #Add axes to figure(s)
    axlist=[self.get_nested(figattr).add_subplot(nrows,ncols,index,**kwargs) for figattr in figlist]
    #Set labels
    for p in ['xlabel','ylabel']:
      if hasattr(axprops,p):
        lbl=getattr(axprops,p)
        if isinstance(lbl,str):
          lblkwd={p:lbl}
        else:
          lblkwd=lbl
        for ax in axlist:
          getattr(ax,'set_'+p)(**lblkwd)
    #Set title
    if hasattr(axprops,'title'):
      lbl=getattr(axprops,'title')
      if isinstance(lbl,str):
        lblkwd={p:lbl}
      else:
        lblkwd=lbl
      for ax in axlist:
        getattr(ax,'set_title')(**lblkwd)
    #Set limits
    for p in ['xlim','ylim']:
      if hasattr(axprops,p):
        vals=getattr(axprops,p)
        for ax in axlist:
          getattr(ax,'set_'+p)(*vals)
    #Store axes list
    self.set_nested(axattr,axlist)
  def iteraxes(self,axlist):
    """Return an iterator over the axes"""
    if axlist is None:
      axlist = [ax for ax in self._all_axattr]
    for axpath in axlist:
      axobj=self.get_nested(axpath)
      #If we got a list of axes, return each of them individually
      if isinstance(axobj,list):
        for ax in axobj:
          yield ax
      else:
        yield axobj
  def axmethod(self,method,axlist=None,outattrs=None,**kwargs):
    """Call the given matplotlib Axes method on the list of axes, with the provided arguments

    Arguments:

      - method = name of Axes method to call, as string
      - axlist = list of nested attribute paths to the axes on which to apply this command, defaults to apply it over all axes
      - outattrs = list of attribute paths to store the results, defaults to None, to not store output
      - \**kwargs = keyword arguments to the method"""
    for idx,ax in enumerate(self.iteraxes(axlist)):
      f=getattr(ax,method)
      if kwargs is None:
        kwd={}
      else:
        kwd=kwargs
      out=f(**kwargs)
      if outattrs is not None:
        self.set_nested(outattrs[idx],out)
    return
  def add_series(self,seriespath,axlist=None,fmt=None,newlabel=None,**kwargs):
    """Add a single series to the listed axes

    Arguments:

      - seriespath = nested path to the series instance to add to the axes
      - axlist = list of nested paths to the Axes instances, defaults to apply it to all axes
      - fmt = matplotlib format string
      - newlabel = label to use for the legend instead of the series-provided one
      - \**kwargs = other keyword arguments for Axes.plot"""
    series=self.get_nested(seriespath)
    for ax in self.iteraxes(axlist):
      series.add_to_axes(ax,fmt,newlabel,**kwargs)
    return
  def add_seriesprops(self,seriesprops):
    """Add a series based on a SeriesProperties instance"""
    #Axes attribute paths
    if hasattr(seriesprops,'axattr'):
      assert not hasattr(seriesprops,'axlist'), "SeriesProperties instance has both axattr (%s) and axlist (%s)"%(str(seriesprops.axattr),str(seriesprops.axlist))
      axlist=[seriesprops.axattr]
    else:
      if hasattr(seriesprops,'axlist'):
        axlist=seriesprops.axlist
      else:
        axlist=None
    #Arguments dictionary
    kwargs=getattr(seriesprops,'kwargs',{})
    for k in seriesprops._direct_attrs:
      if hasattr(seriesprops,k):
        kwargs[k]=getattr(seriesprops,k)
    kwargs['seriespath']=seriesprops.seriesattr
    #Add the series
    self.add_series(axlist=axlist,**kwargs)
  def add_multi_series(self,serieslist,axlist=None,fmtlist=None,labellist=None):
    """Add multiple series to the listed axes

    Arguments:

      - serieslist = list of nested paths to series instances to add to the axes
      - axlist = list of nested paths to Axes instances, defaults to ["ax"]
      - fmtlist = list of format specifier strings, defaults to ``None`` for each series
      - labellist = list of alternative labels for each """
    if fmtlist is None:
      fmtlist=[None]*len(serieslist)
    assert len(serieslist)==len(fmtlist), "Got %d formats for %d series"%(len(fmtlist),len(serieslist))
    if labellist is None:
      labellist=[None]*len(serieslist)
    assert len(serieslist)==len(labellist), "Got %d labels for %d series"%(len(labellist),len(serieslist))
    for i,seriespath in enumerate(serieslist):
      self.add_series(seriespath,axlist,fmtlist[i],labellist[i],**kwargs)
    return
  def hline(self,valpath,axlist=None,**kwargs):
    """Add a horizontal line to the axes

    Arguments:

      - valpath = nested path to the y-value for the horizontal line
      - \**kwargs = keyword arguments for ax.axhline"""
    yval=self.get_nested(valpath)
    for ax in self.iteraxes(axlist):
      ax.axhline(yval,**kwargs)
    return
  def vline(self,valpath,axlist=None,**kwargs):
    """Add a horizontal line to the axes

    Arguments:

      - valpath = nested path to the x-value for the vertical line
      - \**kwargs = keyword arguments for ax.axvline"""
    xval=self.get_nested(valpath)
    for ax in self.iteraxes(axlist):
      ax.axvline(xval,**kwargs)
    return
  def equalityline(self,serpath,axlist=None,label=None,metadata=None):
    """Create a series for a 1:1 line on the requested axes, based on the axis limits

    Because it's based on the axis limits, you should set those before calling this.

    Arguments:

      - serpath = nested attribute path to store the series at
      - axlist = list of nested paths to Axes instances, defaults to ["ax"]
      - label = label to use for the series
      - metadata = metadata to use for the series
    
    I don't ususally use this because it's harder to get the other data series on top of it."""
    for ax in self.iteraxes(axlist):
      xmin,xmax=ax.get_xlim()
      ymin,ymax=ax.get_ylim()
      vals=(min(xmin,ymin), max(xmax,ymax))
      ser=PlotSeries(xvals=vals,yvals=vals,label=label,metadata=metadata)
      self.set_nested(serpath,ser)
    return
  def create_linspace(self,outattr,*args,**kwargs):
    """Call numpy linspace

    This can be useful for defining contour levels.

    Arguments:

      - outattr = nested path to the output data series
      - *args and **kwargs passed to numpy.linspace"""
    ls=np.linspace(*args,**kwargs)
    self.set_nested(outattr,ls)
    return
  def draw_contours(self,axpath,dfpath,xcol,ycol,zcol,levels,
                    vmin=None,vmax=None,cmap="viridis",extend="neither",outattr=None,**kwargs):
    """Create a contour plot.

    Note that for data series to be shown on top of this contour,
    the ``draw_contours`` entry should be in ``basecommands`` instead of ``plotcommands``.

    Arguments:

      - axpath = attribute path to the axes to use for the contour plot 
      - dfpath = attribue path to dataframe containing the data for contouring
      - xcol = name of dataframe column containing x-values (or Stored)
      - ycol = name of dataframe column containing y-values (or Stored)
      - zcol = name of dataframe column containing z-values (or Stored)
      - levels = contour levels (or Stored, for example to result from create_linspace)
      - vmin, vmax = optional (separately) values controling the color scales
      - cmap = color map name, as string
      - extend = which end of the color bar to extend, if either
      - outattr = attribute path for storing the contour data (useful for creating a bar scale later)
      - **kwargs = additional keyword arguments for tricontourf
    """
    #Get the data
    ax=self.get_nested(axpath)
    df=self.get_nested(dfpath)
    xcol=self.get_stored(xcol)
    ycol=self.get_stored(ycol)
    zcol=self.get_stored(zcol)
    levels=self.get_stored(levels)
    #Generate the plot
    cntr = ax.tricontourf(df[xcol], df[ycol], df[zcol], levels, cmap=cmap, vmin=vmin, vmax=vmax,extend=extend,**kwargs)
    ax.axis('equal')
    #Store contour data, if requested
    if outattr is not None:
      self.set_nested(outattr,cntr)
    #Done
    return
  def draw_colorbar(self,axpath,cntrpath,label=None,**kwargs):
    """Draw the color scale bar for a previous contour plot

    Arguments:

      - axpath = attribute path for the axes on which to draw the color bar
      - cntrpath = attribute path for the contour data (from draw_contours, for example)
      - label = color bar label text, as string
      - **kwargs = additional arguments for fig.colorbar"""
    ax=self.get_nested(axpath)
    cntr=self.get_nested(cntrpath)
    label=self.get_stored(label)
    #Draw the color bar
    cbar=fig.colorbar(cntr, ax=ax, **kwargs)
    #Set the label
    if label is not None:
      cbar.set_label(label)
    #Done
    return

#Register for loading from yaml
yaml_manager.register_classes([FigureProperties, AxesProperties, SeriesProperties, FigureRequest])
