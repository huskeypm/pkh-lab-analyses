set size 1,1
#set terminal postscript eps size 4,2.5 solid color enhanced lw 3.0 "Times-Roman" 22 
set terminal postscript eps size 4,3 solid color enhanced lw 3.0 "Times-Roman" 22 

set termoption enhanced
set encoding iso_8859_1
set key top left samplen 1 font ",35"
set macros
POS = "at graph 0.8,0.9 font ',40'"
POL = "at graph 0.5,0.9 font ',35'"
set palette defined (-1 "cyan", 0 "white", 1 "red")

set output "corr-matrix.ps"
set multiplot 

set style fill transparent solid 0.3 noborder

set yrange [0:975]
set xrange[0:972]
set xtics("36" 0, "198" 162, "360" 324, "522" 486, "684" 648, "846" 810, "1007" 972)
set ytics("36" 0, "198" 162, "360" 324, "522" 486, "684" 648, "846" 810, "1007" 972)
#set arrow from 324,0 to 324,972 lc rgb "black" lw 10
set arrow from 324, graph 0 to 324, graph 1 nohead

unset colorbox
#set lmargin 10
#set rmargin 5
################

set origin 0.0,1.0
set size 1,1
#set format x " "
set ylabel "Residue #"
#set xlabel "Residue #"
set title '3 ATP-3 Mg^{2+}' font ",30"
plot '3atp-3mg.dat' matrix with image notitle

set origin 0.9,1.0
set size 1,1
set ylabel " "
#unset ytics
set format y " "
set title '2 ATP-3 Mg^{2+}'
plot '2atp-3mg.dat' matrix with image notitle

set colorbox
#unset lmargin
#unset rmargin
set size 1.1,1.0
set origin 1.8,1.0
set title '2 ATP-2 Mg^{2+}'
plot '2atp-2mg.dat' matrix with image notitle

reset

set origin 0.4,0.0
set arrow from 324,0 to 324,972 nohead lw 3
set size 1,1
#set format x "%g"
set xrange [0:975]
set xtics("36" 0, "198" 162, "360" 324, "522" 486, "684" 648, "846" 810, "1007" 972)
set ytics("36" 0, "198" 162, "360" 324, "522" 486, "684" 648, "846" 810, "1007" 972)
set yrange [0:975]
set xlabel "Residue #" font ",25" 
set ylabel "Residue #" font ",25" 
set title '1 ATP-3 Mg^{2+}' font ",30"
set palette defined (-1 "cyan", 0 "white", 1 "red")
unset colorbox
#set lmargin 10
#set rmargin 5

plot '1atp-3mg.dat' matrix with image notitle

set origin 1.3,0.0
set size 1.1,1
set ylabel " "
#unset ytics
set title '1 ATP-1 Mg^{2+}'
set colorbox
#set lmargin 5
plot '1atp-1mg.dat' matrix with image notitle

unset multiplot
