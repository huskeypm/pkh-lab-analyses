set files { "1atp-1" "1atp-2" "1atp-3" "1atp.1mg-1" "1atp.1mg-2" "1atp.1mg-3" "2atp-1"  "2atp-2"  "2atp-3" "2atp.2mg-1" "2atp.2mg-2"  "2atp.2mg-3"  "3atp-1"  "3atp-2"  "3atp-3" }

set k 0

foreach f $files {
	mol new ../new-trajectories/$f.prmtop
	mol addfile ../wrap-trajectories/$f.pdb
	mol addfile /data/jalenciks/P2X4.mg/wrap-trajectories/$f.dcd type dcd first 0 last -1 step 5 filebonds 1 autobonds 1 waitfor all

	atomselect macro c1 {protein and resid 316 and name CA}
	atomselect macro c2 {protein and resid 640 and name CA}
	atomselect macro c3 {protein and resid 964 and name CA}

	set all [atomselect $k "all"]
	
        set C1 [atomselect $k "c1"]
        set C2 [atomselect $k "c2"]
        set C3 [atomselect $k "c3"]
        
set dat [open $f.txt w]
set num_steps [molinfo $k get numframes]

for {set frame 0} {$frame < $num_steps} {incr frame} {
	$all frame $frame

	$C1 frame $frame
	$C2 frame $frame
	$C3 frame $frame
   
	set dist1 [veclength [vecsub [measure center $C1] [measure center $C2]]]
	set dist2 [veclength [vecsub [measure center $C1] [measure center $C3]]]
	set dist3 [veclength [vecsub [measure center $C2] [measure center $C3]]]

        puts $dat "$frame $dist1 $dist2 $dist3"
} 
set k [ expr $k+1 ]                                                  
}
						             
quit
