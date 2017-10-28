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