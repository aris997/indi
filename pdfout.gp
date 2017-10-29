set terminal postscript enhanced color
set output "Cdt1001_pp11.ps"
set grid
set key
set logscale xy
unset border
set multiplot
set title "Errore di integrazione al variare di dt"
set xlabel 'dt'
set ylabel 'C/Co - 1'
plot 'costantedt01.dat' w l lt rgb 'red'
replot 'costantedt10.dat' w lp lt rgb 'blue'