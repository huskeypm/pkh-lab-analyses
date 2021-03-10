#!/bin/bash


# periodic condition treatment (to move molecules back to box, preventing broken bonds, etc) 

# 1) making sure the molecule is not showing as broken

echo 0 0 | gmx trjconv -f prodrun.xtc -s prodrun.tpr -pbc whole -o tprod1.xtc
echo 0 0 | gmx trjconv -f prodrun2.part0002.xtc -s prodrun.tpr -pbc whole -o tprod2.xtc
echo 0 0 | gmx trjconv -f prodrun2.part0003.xtc -s prodrun.tpr -pbc whole -o tprod3.xtc


# 2) puting the molecule inside the box
echo 0 0 | gmx trjconv -f tprod1.xtc -s step3_input.gro -pbc nojump -o nprod1.xtc
echo 0 0 | gmx trjconv -f tprod2.xtc -s step3_input.gro -pbc nojump -o nprod2.xtc
echo 0 0 | gmx trjconv -f tprod3.xtc -s step3_input.gro -pbc nojump -o nprod3.xtc



rm -r t*xtc

# when viewing in VMD, we will use the step3_input.gro and the new1.xtc

