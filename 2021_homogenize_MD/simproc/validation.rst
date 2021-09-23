Validation
################################################################################


Unit Testing
============

The code in the ``simproc`` package contains a number of doctests to test its components.
These can all be run from the command line with the command:

``python -m simproc --validate``


Integration Testing
===================

The ``validation`` folder contains request files to test a wide variety of the request types
provided by the code.
These tests can be run from the command line with a command such as"

``python -m simproc validation/all_validation.yaml``