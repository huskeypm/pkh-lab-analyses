for j in 2atp 2atp.2mg 1atp 1atp.1mg

do
cpptraj -p ../new-trajectories/$j-1.prmtop -i $j.in
done
