
For usage of singularity, see:
https://sylabs.io/singularity/

Two recipes are provided here:

- `minimal`: contains only the basic software needed to run most simulations
- `complete`: will produce a larger image containing additional useful software for post-processing steps and generation of the documentation

The recipes have been set up so that the resulting singularity images will provide a shell when started with `singularity run <image name>`.

The simulations were originally conducted using fenics 2019.1.
The fenics team no longer provides the necessary data to install fenics 2019.1 in a singularity image.
Instead, the recipes provided here now use fenics 2019.2 instead.
