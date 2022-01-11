set files { "1atp-1" "1atp-2" "1atp-3" "1atp.1mg-1" "1atp.1mg-2" "1atp.1mg-3" "2atp-1"  "2atp-2"  "2atp-3" "2atp.2mg-1" "2atp.2mg-2"  "2atp.2mg-3"  "3atp-1"  "3atp-2"  "3atp-3" }

set k 0

foreach f $files {
        mol new ../new-trajectories/$f.prmtop
        mol addfile ../wrap-trajectories/$f.pdb
        mol addfile /data/jalenciks/P2X4.mg/wrap-trajectories/$f.dcd type dcd first 0 last -1 step 5 filebonds 1 autobonds 1 waitfor all

	set all [atomselect $k "all"]

	set pb [atomselect $k "protein and backbone"]
        set rpb [atomselect $k "protein and backbone" frame 0]
	
	set tmbb [atomselect $k "protein and ((resid 2 to 17) or (resid 300 to 323) or (resid 324 to 342) or (resid 622 to 648) or (resid 649 to 668) or (resid 947 to 972)) and backbone"]
       set rtmbb [atomselect $k "protein and ((resid 2 to 17) or (resid 300 to 323) or (resid 324 to 342) or (resid 622 to 648) or (resid 649 to 668) or (resid 947 to 972)) and backbone" frame 0]

	set p(1) [atomselect $k "pfrag 2 and backbone"]
	set rp(1) [atomselect $k "pfrag 2 and backbone" frame 0]

	set p(2) [atomselect $k "pfrag 1 and backbone"]
	set rp(2) [atomselect $k "pfrag 1 and backbone" frame 0]

        set p(3) [atomselect $k "protein and (resid 649 to 972) and backbone"]
        set rp(3) [atomselect $k "protein and (resid 649 to 972) and backbone" frame 0]

        set pta(1) [atomselect $k "pfrag 2 and resid 2 to 17 and backbone"]
        set pta(2) [atomselect $k "pfrag 1 and resid 324 to 342 and backbone"]
        set pta(3) [atomselect $k "protein and resid 649 to 668 and backbone"]

        set rpta(1) [atomselect $k "pfrag 2 and resid 2 to 17 and backbone" frame 0]
        set rpta(2) [atomselect $k "pfrag 1 and resid 324 to 342 and backbone" frame 0]
        set rpta(3) [atomselect $k "protein and resid 649 to 668 and backbone" frame 0]

        set ptb(1) [atomselect $k "pfrag 2 and resid 300 to 323 and backbone"]
        set ptb(2) [atomselect $k "pfrag 1 and resid 622 to 648 and backbone"]
        set ptb(3) [atomselect $k "protein and resid 947 to 972 and backbone"]

        set rptb(1) [atomselect $k "pfrag 2 and resid 300 to 323 and backbone" frame 0]
        set rptb(2) [atomselect $k "pfrag 1 and resid 622 to 648 and backbone" frame 0]
        set rptb(3) [atomselect $k "protein and resid 947 to 972 and backbone" frame 0]

set dat [open $f-wrap.txt w]
set num_steps [molinfo $k get numframes]

for {set frame 0} {$frame < $num_steps} {incr frame} {
    	$all frame $frame
    	$pb frame $frame
   	$tmbb frame $frame
	
    for {set i 1} {$i < 4} {incr i} {
        $p($i) frame $frame
        $pta($i) frame $frame
        $ptb($i) frame $frame      }

	set trans [measure fit $pb $rpb]
   	$all move $trans
    	set rmsdbb [measure rmsd $pb $rpb]

	set trans [measure fit $tmbb $rtmbb]
    	$all move $trans
    	set rmsdtmbb [measure rmsd $tmbb $rtmbb]

    for {set i 1} {$i < 4} {incr i} {
          set trans [measure fit $p($i) $rp($i)]
          $all move $trans
          set rmsdp($i) [measure rmsd $p($i) $rp($i)]
                                      }

   for {set i 1} {$i < 4} {incr i} {
          set trans [measure fit $pta($i) $rpta($i)]
          $all move $trans
          set rmsdpta($i) [measure rmsd $pta($i) $rpta($i)]
                                      }

   for {set i 1} {$i < 4} {incr i} {
         set trans [measure fit $ptb($i) $rptb($i)]
         $all move $trans
         set rmsdptb($i) [measure rmsd $ptb($i) $rptb($i)]
                                      }

   for {set i 1} {$i < 4} {incr i} {
          set trans [measure fit $pb $rpb]
          $all move $trans
          set rmsdptaa($i) [measure rmsd $pta($i) $rpta($i)]
                                      }

   for {set i 1} {$i < 4} {incr i} {
         set trans [measure fit $pb $rpb]
         $all move $trans
         set rmsdptbb($i) [measure rmsd $ptb($i) $rptb($i)]
                                      }

        for {set i 1} {$i < 4} {incr i} {
        set alpha($i) [expr 57.2958*(acos([vecdot [lindex [measure inertia $pta($i)] 1 2] [lindex [measure inertia $ptb($i)] 1 2]]))]
        if {$alpha($i) > 90} {set alpha($i) [expr 180-$alpha($i)]}
                                        }

         puts $dat "$frame $rmsdbb $rmsdtmbb $rmsdp(1) $rmsdp(2) $rmsdp(3) $rmsdpta(1) $rmsdpta(2) $rmsdpta(3) $rmsdptb(1) $rmsdptb(2) $rmsdptb(3) $alpha(1) $alpha(2) $alpha(3) $rmsdptaa(1) $rmsdptaa(2) $rmsdptaa(3) $rmsdptbb(1) $rmsdptbb(2) $rmsdptbb(3)"
}  
set k [ expr $k+1 ]                                                  
}
						             
quit
