#!/usr/bin/gnuplot -c

set terminal png 

set ytics nomirror
set y2tics nomirror

set title 'Lx vs. steps'

set output ARG1.'.png'


plot ARG1 u 1:2 w lp ti 'Lx'
