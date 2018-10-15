#!/usr/bin/gnuplot -c

set terminal png 


set format x  "10^{%L}"
set logscale x

set xtics nomirror
set ytics nomirror

set x2tics nomirror
set y2tics nomirror


set title 'Fitting CG parameters for eps and r'

set output 'cg_fitting.png'

set xlabel  'energy'
set x2label 'distance'

set ylabel  'badness'
set y2label 'badness'


plot\
 'ave_fit_e1' u 1:2 w lp lw 2 lc '#0000FF'   ti 'eps 1 (energy)',\
 'ave_fit_e2' u 1:2 w lp lw 2 lc '#0000AA'   ti 'eps 2 (energy)',\
 'ave_fit_s1' u 1:2 w lp lw 2 lc '#00FF00'   ti 'sig 1 (r)' axis x2y1,\
 'ave_fit_s2' u 1:2 w lp lw 2 lc '#00AA00'   ti 'sig 2 (r)' axis x2y1
