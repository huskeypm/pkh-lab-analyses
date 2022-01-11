set files { "1atp-1" "1atp-2" "1atp-3" "1atp.1mg-1" "1atp.1mg-2" "1atp.1mg-3" "2atp-1"  "2atp-2"  "2atp-3" "2atp.2mg-1" "2atp.2mg-2"  "2atp.2mg-3"  "3atp-1"  "3atp-2"  "3atp-3" }

set k 0

foreach f $files {
	mol new ../new-trajectories/$f.prmtop
        mol addfile ../wrap-trajectories/$f.pdb
	mol addfile /data/jalenciks/P2X4.mg/wrap-trajectories/$f.dcd type dcd first 0 last -1 step 5 filebonds 1 autobonds 1 waitfor all

	set all [atomselect $k "all"]
	set PRO [atomselect $k "protein and (resid 24 or resid 348 or resid 672) and name CA"]
    	set rPRO [atomselect $k "protein and (resid 24 or resid 348 or resid 672) and name CA" frame 0]
	
	set PRO2 [atomselect $k "protein and (resid 316 or resid 640 or resid 964) and name CA"]
    	set rPRO2 [atomselect $k "protein and (resid 316 or resid 640 or resid 964) and name CA" frame 0]

        set PRO3 [atomselect $k "protein and (resid 312 or resid 636 or resid 960) and name CA"]
        set rPRO3 [atomselect $k "protein and (resid 312 or resid 636 or resid 960) and name CA" frame 0]

set dat [open $f.txt w]
set num_steps [molinfo $k get numframes]

for {set frame 0} {$frame < $num_steps} {incr frame} {
    $all frame $frame
     
    $PRO frame $frame
    $PRO2 frame $frame
    $PRO3 frame $frame
     
    set trans [measure fit $PRO $rPRO]
    $all move $trans
    set prgyr [measure rgyr $PRO]

    set trans [measure fit $PRO2 $rPRO2]
    $all move $trans
    set prgyr2 [measure rgyr $PRO2]

    set trans [measure fit $PRO3 $rPRO3]
    $all move $trans
    set prgyr3 [measure rgyr $PRO3]

    puts $dat "$frame $prgyr $prgyr2 $prgyr3"
} 
set k [ expr $k+1 ]                                                  
}
						             
quit
