set files { "1atp-1" "1atp-2" "1atp-3" "1atp.1mg-1" "1atp.1mg-2" "1atp.1mg-3" "2atp-1"  "2atp-2"  "2atp-3" "2atp.2mg-1" "2atp.2mg-2"  "2atp.2mg-3"  "3atp-1"  "3atp-2"  "3atp-3" }

set k 0

foreach f $files {
        mol new ../../new-trajectories/$f.prmtop
        mol addfile ../../wrap-trajectories/$f.pdb
        mol addfile /data/jalenciks/P2X4.mg/wrap-trajectories/$f.dcd type dcd first 0 last -1 step 5 filebonds 1 autobonds 1 waitfor all

	set all [atomselect $k "all"]
	set potp20 [atomselect $k "resname MG and (within 6 of pfrag 2) and (within 6 of pfrag 0) and (within 10 of (resid 460 or resid 136 or resid 784))"]
	set potp21 [atomselect $k "resname MG and (within 6 of pfrag 2) and (within 6 of pfrag 1) and (within 10 of (resid 460 or resid 136 or resid 784))"]
	set potp10 [atomselect $k "resname MG and (within 6 of pfrag 1) and (within 6 of pfrag 0) and (within 10 of (resid 460 or resid 136 or resid 784))"]

set dat [open $f.txt w]
set num_steps [molinfo $k get numframes]

for {set frame 0} {$frame < $num_steps} {incr frame} {
    $all frame $frame
    $potp20 update
    $potp21 update
    $potp10 update

    $potp20 frame $frame
    $potp21 frame $frame
    $potp10 frame $frame

    set npotp20 [llength [$potp20 get z]]
    set npotp21 [llength [$potp21 get z]]
    set npotp10 [llength [$potp10 get z]]

    puts $dat "$frame $npotp20 $npotp21 $npotp10"
         				       	}
        set k [ expr $k+1 ]
		}
quit
