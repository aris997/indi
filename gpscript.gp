set term x11 persist
set autoscale
a = -2.21082562376599
plot 'costante.dat' u 1 w p
replot a
