set size 1.0,1.0
set terminal postscript eps size 8,2.5 solid color enhanced lw 2.0 "Times-Roman" 22

set termoption enhanced
set encoding iso_8859_1
set key bottom right font ',25' maxrows 3 samplen 1.5

set macros
POS = "at graph 0.01,0.9 font ',35'"
POL = "at graph 0.35,0.9 font ',50'"
set output "corr-vs-res.ps"
set multiplot 

set yrange [-1:1]
set xrange [36:1007]
set ytics -1,0.2,1
set xlabel " " font ",28"
set format x "   "
set xtics format "  "
set ylabel " " font ",28"

set origin 0.0,2.4
set label 2 'D145(MBS2)' @POL
set arrow from 359,-1 to 359,1 nohead lw 3 lt 0
set arrow from 683,-1 to 683,1 nohead lw 3 lt 0
plot "<cat 3atp-3mg.dat | awk 'NR==110{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'orange' t "3ATP-3 Mg^{2+}", "<cat 2atp-2mg.dat | awk 'NR==110{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'red' t "2ATP-2 Mg^{2+}", "<cat 1atp-1mg.dat | awk 'NR==110{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'blue' t "1ATP-1 Mg^{2+}"

unset key
set origin 0.0,1.6
set label 2 'L351(TM Pore)' @POL
plot "<cat 3atp-3mg.dat | awk 'NR==316{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'orange' t "3ATP-3 MG", "<cat 1atp-1mg.dat | awk 'NR==316{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'blue' t "1ATP-1 Mg^{2+}", "<cat 2atp-2mg.dat | awk 'NR==316{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'red' t "2ATP-2 Mg^{2+}"

set origin 0.0,0.8
set label 2 'K70(MBS2)' @POL
set xlabel "Residue #" font ",45"
set ylabel "Correlation" font ",45" offset 0,13
set xtics("36" 36, "198" 198, "360" 360, "522" 522, "684" 684, "846" 846, "1007" 1007)
plot "<cat 3atp-3mg.dat | awk 'NR==64{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'orange' t "3ATP-3 MG", "<cat 1atp-1mg.dat | awk 'NR==64{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'blue' t "1ATP-1 Mg^{2+}", "<cat 2atp-2mg.dat | awk 'NR==64{print $0}' | xargs -n1" u ($0+35):1 w l lw 2 lc rgb 'red' t "2ATP-2 Mg^{2+}"
###############

unset multiplot
