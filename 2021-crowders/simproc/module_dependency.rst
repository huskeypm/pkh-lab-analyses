
Module Dependency Graph
################################################################################

This is only a listing of the dependencies of these modules on each other,
not on the standard library or external libraries.

Keeping track of this is important in order to help prevent circular dependencies.

- package ``requesthandler``

  - ``filepath``: None
  - ``timing``: None
  - ``yaml_manager``: None
  - ``pickle_manager``: None
  - ``nested``: ``yaml_manager``
  - ``locators``: ``filepath``, ``yaml_manager``
  - ``schema``: ``filepath``, ``yaml_manager``, ``locators``, ``nested``
  - ``logging``: ``filepath``, ``timing``, ``yaml_manager``, ``locators``, ``schema``
  - ``request``: ``filepath``, ``yaml_manager``, ``schema``, ``nested``, ``logging``
  - ``locator_requests``: ``yaml_manager``, ``locators``, ``request``, ``logging``
  - ``requestfile``: ``filepath``, ``yaml_manager``, ``locators``, ``request``
  - ``simultaneous``: ``filepath``, ``yaml_manager``, ``request``, ``locators``, ``logging``
  - ``mpi_run``: ``filepath``, ``request``, ``yaml_manager``, ``simultaneous``, ``locators``, ``logging``
  - ``shell``: ``request``, ``yaml_manager``, ``logging``
  - ``debug``: ``request``, ``yaml_manager``, ``shell``, ``logging``
  - ``comparison``: ``request``, ``yaml_manager``, ``logging``
  - ``cleanup``: ``locators``, ``request``, ``yaml_manager``, ``logging``
  - ``customization``: ``filepath``, ``locators``, ``request``, ``yaml_manager``
  - ``commandseq``: ``yaml_manager``, ``pickle_manager``, ``customization``, ``logging``
  - ``templates``: ``yaml_manager``, ``customization``, ``logging``, ``nested``
  - ``generate``: ``yaml_manager``, ``request``, ``nested``, ``customization``
  - ``joblist``: ``yaml_manager``, ``locators``, ``commandseq``, ``logging``
  - ``__init__``: ``filepath``, ``locators``, ``requestfile``, ``customization``, ``shell``, ``templates``, ``cleanup``, ``comparison, ``debug``
  - ``cmdline``: ``*`` (meaning everything listed in ``__init__``)
  - ``__main__``: ``cmdline``

- package ``meshgen``

  - ``buildgeom``: ``requesthandler``
  - ``gmsh_runner``: ``requesthandler``
  - ``dconv_runner``: ``requesthandler``
  - ``hdf5_conv``: ``requesthandler``, ``dconv_runner``
  - ``onestep``: ``requesthandler``, ``gmsh_runner``, ``dconv_runner``, ``hdf5_conv``

- package ``simulation``

  - ``unitsystem``: None
  - ``meshinfo``: ``requesthandler``
  - ``equationbuilder``: None
  - ``plotseries``: None
  - ``simrequest``: ``requesthandler``, ``meshinfo``, ``equationbuilder``, ``plotseries``
  - ``fickian_homog``: ``requesthandler``, ``meshinfo``, ``equationbuilder``, ``simrequest``, ``common_methods``
  - ``ficks_law``:  ``requesthandler``, ``meshinfo``, ``equationbuilder``, ``simrequest``, ``common_methods``

- package ``postproc``

  - ``collection``: ``requesthandler``
  - ``plotting``: ``requesthandler``, ``plotseries``
