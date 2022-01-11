set size 1,1
set terminal postscript eps size 6,3 solid color enhanced lw 3.0 "Times-Roman" 22 

set termoption enhanced
set encoding iso_8859_1
set key top center maxrows 2 font ",25" samplen 2
set macros
POS = "at graph 0.1,0.9 font ',32'"
POL = "at graph 0.8,0.9 font ',40'"

set output "pot-mbs1.ps"
set multiplot 

set style fill transparent solid 0.3 noborder

set xlabel "Time (ns)" font ",25"
set yrange [0:4]
set xrange [-10:1000]
set ytics 0,1,4

################

set origin 0.0,0.0
set ylabel "#K^{+} Ions in center" font ",28"
set xlabel "Time (ns)" font ",25"
set xtics("0" 0, "200" 200, "400" 400, "600" 600, "800" 800, "1000" 1000)

plot '1atp-mean-sd-potct.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "blue" notitle,\
'2atp-mean-sd-potct.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "red" notitle,\
'3atp-mean-sd-potct.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "orange" notitle,\
'2atp.2mg-mean-sd-potct.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "purple" notitle,\
'1atp.1mg-mean-sd-potct.txt' u 1:($5-$6):($5+$6) with filledcurves lc rgb "green" notitle,\
'1atp-mean-sd-potct.txt' u 1:5 w l smooth cspline lw 2 lc rgb "blue" t "1 ATP-3 Mg^{2+}",\
'2atp-mean-sd-potct.txt' u 1:5 w l smooth cspline lw 2 lc rgb "red" t "2 ATP-3 Mg^{2+}",\
'3atp-mean-sd-potct.txt' u 1:5 w l smooth cspline lw 2 lc rgb "orange" t "3 ATP-3 Mg^{2+}",\
'2atp.2mg-mean-sd-potct.txt' u 1:5 w l smooth cspline lw 2 lc rgb "purple" t "2 ATP-2 Mg^{2+}",\
'1atp.1mg-mean-sd-potct.txt' u 1:5 w l smooth cspline lw 2 lc rgb "green" t "1 ATP-1 Mg^{2+}"
#########################

unset multiplot
