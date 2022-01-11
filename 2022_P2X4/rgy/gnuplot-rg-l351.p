set size 1.0,1.0
set terminal postscript eps size 5,4 solid color enhanced lw 3.0 "Times-Roman" 22

set termoption enhanced
set encoding iso_8859_1

set macros
POS = "at graph 0.05,0.9 font ',32'"
POL = "at graph 0.2,0.9 font ',32'"
set output "rg-l351.ps"
set multiplot

set yrange [0:0.2]
set xrange [3:8]
set bmargin 8
set lmargin 15
set ylabel 'Probability' font ",35"
set style fill transparent solid 0.1
set xlabel 'Radius of gyration ({\305})' font ",35"
set key right font ",30" samplen 2 maxrows 5
width = 0.1
bin(x) = width*floor(x/width)

plot "<cat 3atp-*.txt" u (bin($3)):(1) smooth fnormal w l title '3 ATP-3 Mg^{2+}' lc rgb 'orange' lw 3,\
"<cat 2atp.2mg-*.txt" u (bin($3)):(1) smooth fnormal  w l title '2 ATP-2 Mg^{2+}' lc rgb 'red' lw 3,\
"<cat 1atp.1mg-*.txt" u (bin($3)):(1) smooth fnormal w l title '1 ATP-1 Mg^{2+}' lc rgb 'blue' lw 3,\

unset multiplot
