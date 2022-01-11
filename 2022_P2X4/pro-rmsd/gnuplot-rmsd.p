#!/bin/bash

set size 1,1
set terminal postscript eps size 5,2.5 solid color enhanced lw 3.0 "Times-Roman" 22 

set termoption enhanced
set encoding iso_8859_1
set key bottom center maxrows 2 samplen 1
set macros
POS = "at graph 0.04,0.9 font ',32'"
POL = "at graph 0.4,0.9 font ',40'"

set output "rmsd.ps"
set multiplot 

set style fill transparent solid 0.3 noborder

set xlabel "Time (ns)" font ",25"
set yrange [0:6]
set xrange [-10:1000]

################
set origin 0.0,2.4
set label 1 'A' @POS
set label 2 'Whole Protein' @POL`
set xlabel " " font ",25"
set format x "   "
set xtics format "  "
set ylabel " " font ",28"
plot "<cat 1atp-mean-sd-total-protein.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "blue" notitle,\
'2atp-mean-sd-total-protein.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "red" notitle,\
"<cat 3atp-mean-sd-total-protein.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "orange" notitle,\
'2atp.2mg-mean-sd-total-protein.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "green" notitle,\
"<cat 1atp.1mg-mean-sd-total-protein.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "purple" notitle,\
"<cat 1atp-mean-sd-total-protein.txt" u 1:5 w l smooth cspline lw 2 lc rgb "blue" t "1 ATP-3 Mg^{2+}",\
"<cat 2atp-mean-sd-total-protein.txt" u 1:5 w l smooth cspline lw 2 lc rgb "red" t "2 ATP-3 Mg^{2+}",\
"<cat 3atp-mean-sd-total-protein.txt" u 1:5 w l smooth cspline lw 2 lc rgb "orange" t "3 ATP-3 Mg^{2+}",\
"<cat 2atp.2mg-mean-sd-total-protein.txt" u 1:5 w l smooth cspline lw 2 lc rgb "green" t "2 ATP-2 Mg^{2+}",\
"<cat 1atp-mean-sd-total-protein.txt" u 1:5 w l smooth cspline lw 2 lc rgb "purple" t "1 ATP-1 Mg^{2+}"


set origin 0.0,1.6
set label 1 'B' @POS 
set label 2 'P1' @POL
plot "<cat 1atp-mean-sd-p1.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "blue" notitle,\
'2atp-mean-sd-p1.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "red" notitle,\
"<cat 3atp-mean-sd-p1.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "orange" notitle,\
'2atp.2mg-mean-sd-p1.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "green" notitle,\
"<cat 1atp.1mg-mean-sd-p1.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "purple" notitle,\
"<cat 1atp-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "blue" notitle,\
"<cat 2atp-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "red" notitle,\
"<cat 3atp-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "orange" notitle,\
"<cat 2atp.2mg-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "green" notitle,\
"<cat 1atp-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "purple" notitle


set origin 0.0,0.8
set label 1 'C' @POS 
set label 2 'P2' @POL
set ylabel "Backbone RMSD(\305)" font ",35" offset 0,5
plot "<cat 1atp-mean-sd-p2.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "blue" notitle,\
'2atp-mean-sd-p2.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "red" notitle,\
"<cat 3atp-mean-sd-p2.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "orange" notitle,\
'2atp.2mg-mean-sd-p2.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "green" notitle,\
"<cat 1atp.1mg-mean-sd-p2.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "purple" notitle,\
"<cat 1atp-mean-sd-p2.txt" u 1:5 w l smooth cspline lw 2 lc rgb "blue" notitle,\
"<cat 2atp-mean-sd-p2.txt" u 1:5 w l smooth cspline lw 2 lc rgb "red" notitle,\
"<cat 3atp-mean-sd-p2.txt" u 1:5 w l smooth cspline lw 2 lc rgb "orange" notitle,\
"<cat 2atp.2mg-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "green" notitle,\
"<cat 1atp-mean-sd-p1.txt" u 1:5 w l smooth cspline lw 2 lc rgb "purple" notitle


set origin 0.0,0.0
set label 1 'D' @POS 
set label 2 'P3' @POL
set ylabel " " font ",28"
set xlabel "Time (ns)" font ",35"
#set xtics("0" 0, "200" 1000, "400" 2000, "600" 3000, "800" 4000, "1000" 5000)
set xtics("0" 0, "200" 200, "400" 400, "600" 600, "800" 800, "1000" 1000)

plot "<cat 1atp-mean-sd-p3.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "blue" notitle,\
'2atp-mean-sd-p3.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "red" notitle,\
"<cat 3atp-mean-sd-p3.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "orange" notitle,\
'2atp.2mg-mean-sd-p3.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "green" notitle,\
"<cat 1atp.1mg-mean-sd-p3.txt" u 1:($5-$6):($5+$6) with filledcurves lc rgb "purple" notitle,\
"<cat 1atp-mean-sd-p3.txt" u 1:5 w l smooth cspline lw 2 lc rgb "blue" notitle,\
"<cat 2atp-mean-sd-p3.txt" u 1:5 w l smooth cspline lw 2 lc rgb "red" notitle,\
"<cat 3atp-mean-sd-p3.txt" u 1:5 w l smooth cspline lw 2 lc rgb "orange" notitle,\
"<cat 2atp.2mg-mean-sd-p3.txt" u 1:5 w l smooth cspline lw 2 lc rgb "green" notitle,\
"<cat 1atp-mean-sd-p3.txt" u 1:5 w l smooth cspline lw 2 lc rgb "purple" notitle

unset multiplot
