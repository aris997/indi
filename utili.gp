splot '< tail -n 40000 rho100.dat' u 2:3:4 w l // con questo si plotta la coda del file

//plotta dalla riga 2000 alla riga 4000
plot 'costante.dat' every ::2000::4000 w l

set term x11 persist
splot 'drast0.dat' every ::15000::20000 u 2:3:4 w l


unset autoscale

set term x11 persist 0
set multiplot
set xrange [0.4:1.5]
set yrange [0.4:2.8]
set zrange [-0.5:2]
splot 'rho75.dat' u 2:3:4 w l
replot 'rho100.dat' u 2:3:4 w l
unset multiplot

set autoscale

set term x11 persist 1
set multiplot
plot 'rho75.dat' u 1:4 w l
replot 'rho100.dat' u 1:4 w l
unset multiplot

set term x11 persist 2
set multiplot
plot 'rho75.dat' u 1:3 w l
replot 'rho100.dat' u 1:3 w l
unset multiplot

set term x11 persist 3
set multiplot
plot 'rho75.dat' u 1:2 w l
replot 'rho100.dat' u 1:2 w l
unset multiplot




set noxtics //fa sparire tutta la griglia verticale e i numeri su essa







set terminal postscript enhanced color
set output "traiettoria_xyp11.ps"
set grid
unset key
set size square
unset border
set title "Traiettoria sul Piano"
set xlabel "x"
set ylabel "y"
plot 'output.dat' u 2:3 w l lt rgb "red"




set terminal postscript enhanced color
set output "C_pp11.ps"
set grid
unset key
set xrange [0:1000]
set format y '%.12l'
set size square
unset border
set title "Grafico C su tempo"
set xlabel 'tempo [s]'
set ylabel 'C'
plot 'cost.dat' u 1:2 w l lt rgb 'red'

set terminal postscript enhanced color
set output "CCo_pp11.ps"
set grid
unset key
set xrange [0:1000]
set format y '%.2l'
set size square
unset border
set title "Grafico C/Co - 1 su tempo"
set xlabel 'tempo [s]'
set ylabel 'C/Co - 1'
plot 'cost.dat' u 1:3 w l lt rgb 'red'
