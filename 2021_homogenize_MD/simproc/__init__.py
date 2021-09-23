
#Site packages
import pkg_resources #part of setuptools

#Packge modules
from . import requesthandler
from . import meshgen
from . import simulation
from . import postproc

#Paths to files containing doctests
tutorial_file=pkg_resources.resource_filename(__name__,'tutorial.rst')

#Modules for doctest
##requesthandler.doctest_modules.append()
#Files for doctest
requesthandler.doctest_files.append(tutorial_file)