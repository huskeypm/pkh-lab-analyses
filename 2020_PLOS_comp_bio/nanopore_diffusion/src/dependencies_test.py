
"""Check that the required software is available."""

from __future__ import print_function, division #Python 2 compatibility
import os
import subprocess
import sys

#check required python modules
import pathlib #standard library in python 3 but not in python 2
import pkg_resources #part of setuptools
from ruamel.yaml import YAML
import jinja2
import jsonschema
import numpy
import scipy
import matplotlib.pyplot as plt
import pandas
import fenics as fem

#If no subprocess.DEVNULL, open it
if hasattr(subprocess,'DEVNULL'):
  devnull=subprocess.DEVNULL
  doclose=False
else:
  devnull=open(os.devnull,'w')
  doclose=True

#fenics version check
target_fenics_versions=[2016, 2017, 2018]
fenics_version_msg_template="This code was written for FEniCS major versions '%s'. Detected major version '%d'."
assert fem.DOLFIN_VERSION_MAJOR in target_fenics_versions, fenics_version_msg_template%(target_fenics_version,fem.DOLFIN_VERSION_MAJOR)

#check that gmsh runs
gmsh_cmd="gmsh --version"
retcode=subprocess.call(gmsh_cmd,stderr=devnull,shell=True)
assert retcode==0, "gmsh failed to execute correctly"
if doclose:
  devnull.close()

#tests passed
print("OK")
