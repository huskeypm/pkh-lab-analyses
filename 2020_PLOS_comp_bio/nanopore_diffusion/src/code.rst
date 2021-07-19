Code Documentation
################################################################################

General Overview
================

This collection of scripts can perform the following steps:

#. Generate meshes for FEM simulation: `Mesh Generation`_
#. Run simulations in ``FEniCS`` and extract selected data: Simulation_
#. Post-process the simulation results to generate graphics: `Post-Processing`_
#. Generate lengthy input files for the steps above containing parametric combinations of input values: `Parameter Generation`_

There are other additional `Miscellaneous Modules`_ as well.

**TODO**: list required software, as in README

**TODO**: describe expected folder structure, as in README (including environment variables)

**TODO**: maybe some other stuff to pull in from README as well

Input Files
-----------

The scripts makes extensive use of ``yaml`` format (https://en.wikipedia.org/wiki/YAML) for input files,
and even some output files as well.
A very brief introduction to the syntax can be found at https://learnxinyminutes.com/docs/yaml/,
and its official site is http://yaml.org/.

Each yaml input file can contain one or more documents.
Each document defines a single instance of a class defined by one of the python modules,
usually (but not always) the same module which accepts the yaml file as input.
The yaml files should contain a comment indicating which class each document is converted to.

The attributes of the class are set by the dictionary values defined in each yaml document.
Some of these attributes will be set to numerical values, strings, or other python built-in data types like lists and dictionaries.
Some of them, however, are further processed into other classes defined by the python module.
The documentation of each class will indicate any attributes which are converted to other classes.

Thus, the key to understanding the input files is to look at the appropriate class definitions in this documentation.

In particular, note that many parameters for a simulation are defined in the simulator module it uses,
not in ``simulator_general.py`` or ``simulator_run.py``.

**TODO** give an example of how to interpret a YAML file from the class documentation.

Command-Line Execution
----------------------

Each of the 4 steps in the input process has a script that can be executed from the command line.
The first argument should be the input yaml file.
Generally, these scripts also accept the argument ``--help`` to provide more information on their usage.

**TODO**: explain other arguments supported by the common module.

Task Automation with DoIt
-------------------------

**TODO** this section, including explanation of ``dodo.py``

.. _`Mesh Generation`:
Mesh Generation
===============

This process creates a ``.geo`` file, which ``gmsh`` then processes into ``.msh`` format,
which is then converted to ``.xml`` usable by ``FEniCS`` with ``dolfin-convert``.
But sometimes we want to run code in parallel,
so we also have to use ``FEniCS`` to convert the mesh to HDF5 format.

buildgeom.py
------------

.. automodule:: buildgeom
   :members:

geom_mk_msh.py
--------------

.. automodule:: geom_mk_msh
   :members:

geom_mk_xml.py
--------------

.. automodule:: geom_mk_xml
   :members:

geom_mk_hdf5.py
---------------

.. automodule:: geom_mk_hdf5
   :members:

.. _Simulation:
Simulations
===========

Each type of simulation that can be performed has its own simulator module.
These modules all use clases defined in simulator_general_.
The correct simulator module can be run by simulator_run_.

Individual simulations can also customize the behavior of the simulator by requesting modules from ``customizations``.

**TODO** explain how customization works

.. _simulator_general:
simulator_general.py
--------------------

.. automodule:: simulator_general
   :members:

.. _simulator_run:
simulator_run.py
----------------

.. automodule:: simulator_run
   :members:

fickian_unhomog.py
------------------

.. automodule:: fickian_unhomog
   :members:

smol_unhomog.py
---------------

.. automodule:: smol_unhomog
   :members:

tdpnp_unhomog.py
----------------

.. automodule:: tdpnp_unhomog
   :members:

.. _`Post-Processing`:
Post-Processing
===============

This step is run by postproc.py_, which can call collect_results_ and plotdata_.

.. _postproc.py:
postproc.py
-----------

.. automodule:: postproc
   :members:

.. _collect_results:
collect_results.py
------------------

.. automodule:: collect_results
   :members:

.. _plotdata:
plotdata.py
-----------

.. automodule:: plotdata
   :members:

.. _`Parameter Generation`:
Parameter Generation
====================

paramgen.py
-----------

.. automodule:: paramgen
   :members:

.. _`Miscellaneous Modules`:
Miscellaneous Modules
=====================

common.py
---------

.. automodule:: common
   :members:

folderstructure.py
------------------

.. automodule:: folderstructure
   :members:

dependencies_test.py
--------------------

.. automodule:: dependencies_test
   :members:

unitsystem.py
-------------

.. automodule:: unitsystem
   :members:

