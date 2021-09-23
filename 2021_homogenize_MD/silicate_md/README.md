
# Folder descriptions

- **common**: files that are re-used in the runs
- **jobs**: storage location for MD jobs (empty until jobs are created)
- **modules**: Python modules used by requests and scripts
- **notebooks**: Jupyter notebooks used for testing and investigation
- **requests**: input data files for generating simulation data files from templates and scripts
- **scripts**: Python scripts used in pre or post-processing of the MD data

Within `jobs/prepared`, `jobs/run`,
the simulations are organized into the following levels:
- set (folders tracked in git)
- job

Each MD analysis can be divided into "stages" labeled with a short string.

# Simulation Files

Within each simulation directory, most file names are prefixed with a three-digit number
to help organize the different stages of the analysis.
The first digit in these numbers is as follows:

- 0 or 1: analysis setup files
- 2, 3, 4, 5: MD simulation input/output files
- 6, 7, 8, 9: Files created during post-processing
