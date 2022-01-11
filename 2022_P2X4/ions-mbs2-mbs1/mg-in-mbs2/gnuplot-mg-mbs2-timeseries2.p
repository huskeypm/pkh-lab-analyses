set size 1,1
set terminal postscript eps size 4,2 solid color enhanced lw 3.0 "Times-Roman" 22

set termoption enhanced
set encoding iso_8859_1
set key right center font ',30' maxrows 3 samplen 1.5

set macros
POS = "at graph 0.05,0.8 font ',32'"
POL = "at graph 0.1,0.6 font ',25'"
set output "mg-mbs2-timeseries2.ps"
set multiplot 

set xrange [-10:1000]
set yrange [-0.1:1.1]
set ytics 0,1,1

################
set xlabel "Time (ns)" font ",30"
set xtics("0" 0, "200" 200, "400" 400, "600" 600, "800" 800, "1000" 1000)
set ylabel "Mg^{2+} in MBS2" font ",30"

set origin 0.0,0.8
set label 2 'P1-P3' @POL
plot "<cat 3atp-1.txt | awk '$2>0{print}'" u 1:2 w l lw 2 lc rgb "orange" t '3 ATP-3 Mg^{2+}',\
'2atp.2mg-1.txt' u 1:2 w l lw 2 lc rgb "red" t '2 ATP-2 Mg^{2+}',\
"<cat 1atp.1mg-1.txt | awk '$2>0{print}'" u 1:($2+0.05) w l lw 2 lc rgb "blue" t '1 ATP-1 Mg^{2+}'

set origin 0.9,0.8
set ylabel " " font ",25"
set format y "   "
set ytics format "  "
set label 2 'P1-P2' @POL
plot "<cat 3atp-1.txt | awk '$3>0{print}'" u 1:3 w l lw 2 lc rgb "orange" notitle,\
"<cat 2atp.2mg-1.txt | awk '$3>0{print}'" u 1:($3+0.05) w l lw 2 lc rgb "red" notitle,\
"<cat 1atp.1mg-1.txt" u 1:3 w l lw 2 lc rgb "blue" notitle

set origin 1.8,0.8
set label 2 'P2-P3' @POL
plot "<cat 3atp-1.txt | awk '{print}'" u 1:4 w l lw 2 lc rgb "orange" notitle,\
"<cat 2atp.2mg-1.txt | awk '$4>0{print}'" u 1:($4+0.05) w l lw 2 lc rgb "red" notitle,\
"<cat 1atp.1mg-1.txt" u 1:4 w l lw 2 lc rgb "blue" notitle
################

unset multiplot
