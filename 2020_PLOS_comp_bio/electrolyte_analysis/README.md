
Data for simulations.

Makes use of nanoscale diffusion code.
This repo is meant to be used as the datafolder for that code.

# Runs

- run001: Run Shashank's original data
- run002: Set nonzero potential on surface of sphere; do simultaneous simulations
- run003: Use multi-species code, with reactive boundary condition; use jinja2 templates; reduce number of cases (not full graph)
- run004: body-centered mesh
- run005: 2D mesh with random circular inclusions
- run006: separate properties of main inclusion from others
- run007: test case for diffusive inclusions
- run008: two reactive boundaries, non-diffusive inclusions
- run009: 3D mesh (body-centered) with two reactive boundaries
- run010: 3D mesh (body-centered) with two reactive boundaries, solving first-order system (not working yet)
- run011: 3D mesh (body-centered) with no reactive boundaries, using Dirichlet conditions instead.
- run012: 3D mesh (body-centered) with sequential solution for reactive boundaries

# Folders
- logs: output from the scripts and other programs when run
- mesh: mesh data files, including templates but not parameter definition files
- model: model template files
- notebooks: jupyter notebooks, used for post-processing
- old: old files no longer used
- params: parameter definition yaml files
- postproc: post-processing outputs
- solutions: simulation output data

# Workflow
At the top level, each run as a script to do the following:
  - generate the input files (.geo files and parameter yaml files)
  - do all the steps of mesh generation (.geo->.msh->.xml->.hdf5)
  - run the simulations
The script puts logs in the "logs" folder.

After the simulations are complete, the post-processing notebook can be
used to collect the data and generate plots.

Two top-level scripts are provided, which take a basename as an argument:
  - do.sh: run the script (described above) for the specified run
  - cleanup.sh: remove all generated output files for the specified run

Prior to running these scripts, set the `DATALOC` and `SRCLOC` environment variables:
`source ./set_env.sh`
You may need to customize this script based on where you store the data and source code.

# To create a new run

1) create its mesh template in `mesh/templates`
2) create its model template in `model/templates`
3) create its input generation script in `scripts`
4) create the shell script to run it in the root folder for the repo.

