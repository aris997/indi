set term x11 persist 0
set multiplot
splot 'rho75.dat' u 2:3:4 w l
replot 'rho100.dat' u 2:3:4 w l
unset multiplot

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
