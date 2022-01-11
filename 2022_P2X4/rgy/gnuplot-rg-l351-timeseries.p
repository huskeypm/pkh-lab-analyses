set size 1.0,1.0
set terminal postscript eps size 8,2.5 solid color enhanced lw 3.0 "Times-Roman" 22

set termoption enhanced
set encoding iso_8859_1
set key right font ',35' maxrows 1 samplen 1.5

set macros
POS = "at graph 0.05,0.9 font ',32'"
POL = "at graph 0.2,0.9 font ',32'"
set output "rg-l351-timeseries.ps"
set multiplot 

set xrange [0:1000]
set yrange [4:9]

################
set xlabel " " font ",28"
set format x "   "
set xtics format "  "

set ylabel " " font ",28"
set ytics 0,1,10

set origin 0.0,4.0
set label 1 'A' @POS
set label 2 '1-ATP' @POL
plot '1atp-1.txt' u 1:3 w l lt 0 lw 2 lc rgb "blue" t 'set-1',\
'1atp-2.txt' u 1:3 w l lc rgb "blue" t 'set-2',\
'1atp-3.txt' u 1:3 w l lt 0 lw 5 lc rgb "blue" t 'set-3'  

set origin 0.0,3.2
set label 1 'B' @POS
set label 2 '1-ATP(1 Mg^{2+})' @POL
plot '1atp.1mg-1.txt' u 1:3 w l lt 0 lw 2 lc rgb "green" t 'set-1',\
'1atp.1mg-2.txt' u 1:3 w l lc rgb "green" t 'set-2',\
'1atp.1mg-3.txt' u 1:3 w l lt 0 lw 5 lc rgb "green" t 'set-3'

set origin 0.0,2.4
set label 1 'C' @POS 
set label 2 '2-ATP' @POL
plot '2atp-1.txt' u 1:3 w l lt 0 lw 2 lc rgb "red" t 'set-1',\
'2atp-2.txt' u 1:3 w l lc rgb "red" t 'set-2',\
'2atp-3.txt' u 1:3 w l lt 0 lw 5 lc rgb "red" t 'set-3'  

set origin 0.0,1.6
set label 1 'D' @POS
set label 2 '2-ATP(2 Mg^{2+})' @POL
plot '2atp.2mg-1.txt' u 1:3 w l lt 0 lw 2 lc rgb "purple" t 'set-1',\
'2atp.2mg-2.txt' u 1:3 w l lc rgb "purple" t 'set-2',\
'2atp.2mg-3.txt' u 1:3 w l lt 0 lw 5 lc rgb "purple" t 'set-3'

set origin 0.0,0.8
set label 1 'E' @POS 
set label 2 '3-ATP' @POL
set xlabel "Time (ns)" font ",40"
set xtics("0" 0, "200" 200, "400" 400, "600" 600, "800" 800, "1000" 1000)
set ylabel "R_{g} of TM residue L351(\305)" font ",40" offset 0,26
plot '3atp-1.txt' u 1:3 w l lt 0 lw 2 lc rgb "orange" t 'set-1',\
'3atp-2.txt' u 1:3 w l lc rgb "orange" t 'set-2',\
'3atp-3.txt' u 1:3 w l lt 0 lw 5 lc rgb "orange" t 'set-3'  

################

unset multiplot
