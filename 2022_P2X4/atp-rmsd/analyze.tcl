set files { "1atp-1" "1atp-2" "1atp-3" "1atp.1mg-1" "1atp.1mg-2" "1atp.1mg-3" "2atp-1"  "2atp-2"  "2atp-3" "2atp.2mg-1" "2atp.2mg-2"  "2atp.2mg-3"  "3atp-1"  "3atp-2"  "3atp-3" }

set k 0

foreach f $files {
        mol new ../new-trajectories/$f.prmtop
        mol addfile ../wrap-trajectories/$f.pdb
        mol addfile /data/jalenciks/P2X4.mg/wrap-trajectories/$f.dcd type dcd first 0 last -1 step 5 filebonds 1 autobonds 1 waitfor all

	set all [atomselect $k "all"]

       set atp(0) [atomselect $k "(resid 975) and noh and resname ATP"]
       set atp(1) [atomselect $k "(resid 976) and noh and resname ATP"]
       set atp(2) [atomselect $k "(resid 977) and noh and resname ATP"]
       set atp(3) [atomselect $k "(resid 978) and noh and resname ATP"]
       set atp(4) [atomselect $k "(resid 974) and noh and resname ATP"]

        set ratp(0) [atomselect $k "(resid 975) and noh and resname ATP" frame 0]
        set ratp(1) [atomselect $k "(resid 976) and noh and resname ATP" frame 0]
        set ratp(2) [atomselect $k "(resid 977) and noh and resname ATP" frame 0]
        set ratp(3) [atomselect $k "(resid 978) and noh and resname ATP" frame 0]
        set ratp(4) [atomselect $k "(resid 974) and noh and resname ATP" frame 0]

set dat [open $f-wrap.txt w]
set num_steps [molinfo $k get numframes]

for {set frame 0} {$frame < $num_steps} {incr frame} {
    	$all frame $frame
	
    for {set i 0} {$i < 5} {incr i} {
        $atp($i) frame $frame      
				    }

   for {set i 4} {$i < 5} {incr i} {
   	set trans [measure fit $atp($i) $ratp($i)]
   	$all move $trans
   	set rmsdatp($i) [measure rmsd $atp($i) $ratp($i)]
                                   }

         puts $dat "$frame $rmsdatp(0) $rmsdatp(1) $rmsdatp(2) $rmsdatp(3) $rmsdatp(4)"
}  
set k [ expr $k+1 ]                                                  
}
						             
quit
