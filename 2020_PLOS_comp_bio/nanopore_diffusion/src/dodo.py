"""Doit file for model runs

Discussion:

This file is intended to automate the process from mesh generation, through running the simulator, through post-processing.
As such, it relies on scripts specific to each of those to generate the necessary tasks.
The runs requested in control.yaml are used."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp

#Site packages
from doit import create_after

#Local
import folderstructure as FS
import common
import paramgen
import buildgeom
import geom_mk_msh
import geom_mk_xml
import geom_mk_hdf5
import simulator_run
import postproc

#Constants
controlfile=osp.join(FS.datafolder,'control.yaml')

#Read in the files for processing
infile_list=common.readyaml(controlfile)


def generic_task_generator(folder,objtype,infile_list=infile_list):
  for infile in infile_list:
    infpath=osp.join(folder,infile)
    if osp.isfile(infpath):
      print("Loading %s from %s."%(objtype.__name__,infpath))
      allobj=objtype.all_from_yaml(infpath)
      for obj in allobj:
        tdo=obj.task_definition
        if hasattr(tdo,'__next__'):
          #Got a generator: yield everything it does
          for td in tdo:
            yield td
        else:
          #Single task definition to yield
          yield tdo

#Parameter generation tasks
def task_paramgen():
  return generic_task_generator(FS.params_paramgen_folder,paramgen.ParameterGenerator)

#Mesh tasks
@create_after('paramgen')
def task_make_geo():
  return generic_task_generator(FS.params_mesh_folder,buildgeom.MeshParameters)

@create_after('paramgen')
def task_make_msh():
  return generic_task_generator(FS.params_mesh_folder,geom_mk_msh.GmshRunner)

@create_after('paramgen')
def task_make_xml():
  return generic_task_generator(FS.params_mesh_folder,geom_mk_xml.DolfinConvertRunner)

@create_after('paramgen')
def task_make_hdf5():
  return generic_task_generator(FS.params_mesh_folder,geom_mk_hdf5.HDF5Converter)

#Simulator tasks
@create_after('paramgen')
def task_solve():
  return generic_task_generator(FS.params_model_folder,simulator_run.ModelParameters)

#Post-processing tasks
@create_after('paramgen')
def task_postproc():
  return generic_task_generator(FS.params_postproc_folder,postproc.PostProcParameters)
