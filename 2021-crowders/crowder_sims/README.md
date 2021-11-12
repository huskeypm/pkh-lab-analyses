
Simulation files for crowder simulations.

# Folders

- logs: output from the scripts and other programs when run
- mesh: mesh data files, including templates but not parameter definition files
- modules: python code to be imported
- postproc: post-processing outputs
- requests: request input files
- scripts: supplemental python scripts
- solutions: simulation output data

# Requests

- cleanup: remove output of
- test01: basic functionality testing with only explicit inclusions
- test02: testing generation of inclusions on a regular lattice
- test03: test a suite of analyses, with a regular lattice of crowders
- test04: test increasing numbers of crowders by reducing separation buffer (_do not use_)
- test05: allow alteration of dirichlet condition on crowders (_crashes, do not use_)
- test06: try MPI for fenics steps
- test07: add loop for crowder potential
- test08: partially reactive inclusions
- test09: same simulations as test07, but using job list instead of template loops
- test10: same simulations as test08, but using job list instead of template loops
- test11: single reactive inclusion, vary its potential and area fraction
- test12: single reactive inclusion, variable boundary condition on pore
- test13: perturb crowder positions and run repeated iterations

# Setup

The files in this repo depend on the `nanopore_diffusion` repo,
which should be added to the python path.
This can be done without creating copies of the `nanopore_diffusion` repo.
From within the `nanopore_diffusion` folder, execute:
`python setup.py develop --user`
This is only needed once.

# Usage

If using the singularity image, start it with
`singularity run fenics_2019.simg`
from the directory containing the image file,
and then change into this directory.

A script is provided here to set the necessary environment variables before running.
`source ./set_env.sh`

To run the requests in a given file:
`python -m simproc <request file>`
or, with doit:
`doit control=<request file>`

Some requests may generate other request files.
The generated request files are not run automatically.
A separate command is required to run them.
For example, to run test03, both of the following steps are required, in order:
`python -m simproc requests/mk_test03.yaml`
`python -m simproc requests/generated/test03.yaml`
There is a convenience script for this at `scripts/runjob.sh`:
`runjob.sh test03` will perform both of the actions above.

To clean up the output from a particular request file, use the `--select` argument:
`python -m simproc requests/cleanup.yaml --select clean_logs test01.yaml`

To run a jupyter notebook under singularity:
`export XDG_RUNTIME_DIR=""`
`jupyter notebook --no-browser`

# Workflow

Mesh Generation:

- a gmsh input file (.geo) is generated from a template (.geo.jinja2).
- gmsh is run to create the .msh file from the .geo file. (Other output files are also created.)
- dolfin-convert is run to create an xml file from the msh file.
- A module in `nanopore_diffusion` is called to convert the xml file to hdf5 format.

FEM Simulations:

- The LPB equation is solved to compute the potential
- The Smoluchowski equation is solved to compute concentration fields
- The data is post-processed