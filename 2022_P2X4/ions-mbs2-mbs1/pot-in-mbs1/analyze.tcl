set files { "1atp-1" "1atp-2" "1atp-3" "1atp.1mg-1" "1atp.1mg-2" "1atp.1mg-3" "2atp-1"  "2atp-2"  "2atp-3" "2atp.2mg-1" "2atp.2mg-2"  "2atp.2mg-3"  "3atp-1"  "3atp-2"  "3atp-3"}

set k 0

foreach f $files {
        mol new ../new-trajectories/$f.prmtop
        mol addfile ../wrap-trajectories/$f.pdb
        mol addfile /data/jalenciks/P2X4.mg/wrap-trajectories/$f.dcd type dcd first 0 last -1 step 5 filebonds 1 autobonds 1 waitfor all

	set all [atomselect $k "all"]
	set potct [atomselect $k "element K and (within 8 of pfrag 0) and (within 8 of pfrag 1) and (within 8 of pfrag 2) and (not within 10 of lipids)"]

set dat [open $f.txt w]
set num_steps [molinfo $k get numframes]

for {set frame 0} {$frame < $num_steps} {incr frame} {
    $all frame $frame
    $potct update
    #$potct frame $frame

    set npotct [llength [$potct get z]]

    puts $dat "$frame $npotct"
         				       	}
        set k [ expr $k+1 ]
		}

quit
