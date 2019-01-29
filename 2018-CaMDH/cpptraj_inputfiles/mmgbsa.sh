#!/bin/bash
#
# This script is used to calcluate the MMGBSA binding free enengy 
# between DH and CaM/CaMBR.
#

# The sampling (frames) of DH, CaM/CaMBR and DH-CaM/CaMBR were extracted
# from the same MD trajectories..


# the temp.binpos is new trajectory file generated 
# by extracting frame every 2 ns from original trajectory
# Thus the MMGBSA calculation will be pefromed at about 500 frames 
# (the cumulative simlulation length is ~1 us for each case)



# 1) extract DH trajectory (which is the ligand), the residue is from :196-213
cat > extract_dh.trajin << EOD
trajin temp.binpos
strip !:196-213
trajout DH.traj
EOD
cpptraj -p x.top -i extract_dh.trajin


# 2) extract CaM/CaMBR trajectory (which is the receptor), the residue is :1-168
cat > extract_camcambr.trajin << EOD
trajin temp.binpos
strip !:1-168
trajout camcambr.traj
EOD
cpptraj -p x.top -i extract_camcambr.trajin


# 3) extract the CaM/CaMBR + DH (which is the complex). the resideu is :1-168 || :196-213
cat > extract_complex.trajin << EOD
trajin temp.binpos
strip !(:1-168,196-213)
trajout complex.traj
EOD
cpptraj -p x.top -i extract_complex.trajin


# 4) now generete the corresponding topology files for DH, camcambr and complex.
cat > temp_trajin << EOD
parm x.top
parmstrip !:196-213
parmwrite out DH.parm7
EOD

cpptraj -i temp_trajin

cat > temp_trajin << EOD
parm x.top
parmstrip !:1-168
parmwrite out camcambr.parm7
EOD

cpptraj -i temp_trajin

cat > temp_trajin << EOD
parm x.top
parmstrip !(:1-168,196-213)
parmwrite out complex.parm7
EOD

cpptraj -i temp_trajin

# Now performe the MMGBSA calculation
cat > mmgbsa.in << EOD
&general
verbose=2, keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=2,
/
EOD


mpirun -np 6 MMPBSA.py.MPI -O -i mmgbsa.in -cp complex.parm7 -rp camcambr.parm7 -lp DH.parm7 -y complex.traj -yr camcambr.traj -yl DH.traj

