
#!/bin/bash


#The high linker flexibility  CaM model in which the entire linker was subject to a short annealing MD in vacuum by Amber to destroy the $\alpha$-helix structure with DSSP predicted secondary structure as coils and bend ("CCCSCCSSSCCCSSS" where "C" and "S" refers to coil and bend, respectively).

# NOTE: dssp does not recognize HIE (amber type)!!! need to change this to HIS 
python martini.py  -f 3cln_linkerCoil.pdb -o 3cln_CG.top -x 3cln_CG.pdb -dssp /usr/bin/dssp -p backbone -ff martini22

# -- for the normal linker flexibility, replace 3cln_linkerCoil.pdb with 3cln.pdb

# do this for the CaMBR
python martini.py  -f cambr.pdb -o cambr_CG.top -x cambr_CG.pdb -dssp /usr/bin/dssp -p backbone -ff martini22

# place the CaMBR about 6nm away from CaM
gmx insert-molecules -f 3cln_CG.pdb -ci cambr_CG.pdb -box 6 6 6 -nmol 1 -radius 0.21 -try 100 -o combine1.gro



# change  the top file to reflect the inclusion of CaMBR

# minimize the side chain in vacuum
gmx grompp -f minimization.mdp -p system.top -c combine1.gro -o vacummMin.tpr -r combine1.gro 
gmx mdrun -deffnm vacummMin -v


# put the molecule at the center of the box
gmx editconf -f combine1.gro -d 2.5 -bt cubic -c -o centeredstr.gro


# add waters
  gmx solvate -cp centeredstr.gro -cs water.gro -box 15 15 15 -radius 0.21 -o solvated.gro


gmx grompp -f minimization.mdp -p system.top -c solvated.gro -o solvtMin.tpr -r solvated.gro 
gmx mdrun -deffnm solvtMin -v


# before adding ions, get the number of waters

nW=`tail -n 1 system.top | awk '{ print $2 }'`


#add Na and get the number of replaced waters

gmx insert-molecules -f solvated.gro -ci Na.gro -nmol 305 -replace 'resname W' -radius 0.21 -try 500 -o withNa.gro 2> t1.log
W1=`grep Replaced t1.log | awk '{ print $2 }'`


# add Cl and get the number of replaced waters
gmx insert-molecules -f withNa.gro -ci Cl.gro -nmol 305 -replace 'resname W' -radius 0.21 -try 500 -o withClNa.gro 2> t2.log

W2=`grep Replaced t2.log | awk '{ print $2 }'`

# set the right number of waters in the top file
cp temp.top system.top
newW=`echo | awk '{ print '$nW' -  '$W1' - '$W2' }'`
sed -i 's/XXX/'$newW'/g' system.top


# minimization
gmx grompp -f minimization.mdp -p system.top -c withClNa.gro -o withClNaMin.tpr -r withClNa.gro 2>t5log
gmx mdrun -deffnm withClNaMin -v 2> t6log


