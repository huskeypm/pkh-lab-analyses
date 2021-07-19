"""Post-processing to generate additional data and graphics"""

#Standard library
import os.path as osp

#Site packages

#Local
import folderstructure as FS
import simulator_general
import plotdata
import collect_results
import common

#Constants

class PostProcParameters(common.ParameterSet):
  """Definition of post-processing requests

  Attributes:

    - modelparamsfile = name of yaml file containing the models to post-process (include .yaml extension)
    - modelnames = sequence of model names that this post-processing definition applies to

      Omit or set to None to apply to all model names for the basename.
      An empty list will result in no post-processing being done.

    - do_collection = boolean, True to collect output infofiles into a Pandas DataFrame (see folderstructure for definition of infofile)
    - collection_exclusions = list of keys in output infofile to be excluded from DataFrame (see folderstructure for definition of infofile)
    - model_plots = sequence of plots to generate for each model

      Each plot defines a plotdata.PlotFigure instance.

    - collection_plots = sequence of plots to generate from the collected DataFrame

      Each plot defines a plotdata.PlotFigure instance.
      If collection_schema is omitted or None, this should be as well, or be an empty sequence.

  Attributes to be created by methods:

    - modellist = list of ModelParametersBase instantaces to postprocess, """
  __slots__=['modelparamsfile','modelnames','do_collection','collection_exclusions','model_plots','collection_plots','modellist']
  _required_attrs=['basename','modelparamsfile']
  #Remember loaded models (parameters, not simulators) so we aren't reading the same files over and over
  loaded_modelfiles={} #{model filename: [model names]}
  loaded_models={} #{model name: ModelParametersBase instance}

  def __init__(self,**kwd):
    #Initialization from base class
    super(PostProcParameters, self).__init__(**kwd)
    #Load the model parameters file if not already loaded
    if not self.modelparamsfile in self.loaded_modelfiles.keys():
      self.loaded_modelfiles[self.modelparamsfile]=[]
      modelparams_fpath=osp.join(FS.params_model_folder,self.modelparamsfile)
      modelparams_gen=simulator_general.ModelParametersBase.all_from_yaml(modelparams_fpath)
      for mp in modelparams_gen:
        #Make sure this modelname is not already used
        if mp.modelname in self.loaded_models.keys():
          for fname,mlist in self.loaded_modelfiles.items():
            if mp.modelname in mlist:
              raise Exception("Duplicate model name: %s in both %s and %s."%(mp.modelname,fname,self.modelparamsfile))
        self.loaded_models[mp.modelname]=mp
        self.loaded_modelfiles[self.modelparamsfile].append(mp.modelname)
    #Get the indicated modelnames
    if getattr(self,'modelnames',None) is not None:
      self.modellist=[v for k,v in self.loaded_models.items() if k in self.modelnames]
    else:
      self.modellist=[self.loaded_models[mname] for mname in self.loaded_modelfiles[self.modelparamsfile]]

  def allsubs(self):
    #Generate model plots, if any
    if getattr(self,'model_plots',None) is not None:
      #Each requested model plot
      for mplotdict in self.model_plots:
        for modelparams in self.modellist:
          mplotdict.update({'basename':self.basename,'modelname':modelparams.modelname})
          mplot = plotdata.ModelPlotFigure(**mplotdict)
          yield mplot
    #Do the collection, if requested
    if getattr(self,'do_collection',False):
      collector=collect_results.ResultsCollector(basename=self.basename,modellist=self.modellist,exclusions=getattr(self,'collection_exclusions',[]))
      yield collector
      #Generate each requested plot from the collected data
      if getattr(self,'collection_plots',None) is not None:
        for cplotdict in self.collection_plots:
          cplotdict.update({'basename':self.basename})
          cplot=plotdata.CollectionPlotFigure(**cplotdict)
          yield cplot

  @property
  def task_definition(self):
    for obj in self.allsubs():
      yield obj.task_definition

  def run(self):
    for obj in self.allsubs():
      obj.run()

#Support command-line arguments
if __name__ == '__main__':
  program_description='Run post-processing commands'
  input_file_description="Path to the file containing PostProcParameters definitions"
  
  common.run_cmd_line(program_description,input_file_description,PostProcParameters)
