set term x11 persist 0
set multiplot
unset autoscale
set xrange [0:400]
set yrange [0:3.2]
plot 'drast030.dat' u 1:2 w l
replot 'drast030.dat' u 1:3 w l
replot 'drast030.dat' u 1:4 w l
unset multiplot

set term x11 persist 1
set multiplot
unset autoscale
set xrange [0:400]
set yrange [0:3.2]
plot 'drast050.dat' u 1:2 w l
replot 'drast050.dat' u 1:3 w l
replot 'drast050.dat' u 1:4 w l
unset multiplot

set term x11 persist 2
set multiplot
unset autoscale
set xrange [0:400]
set yrange [0:3.2]
plot 'drast070.dat' u 1:2 w l
replot 'drast070.dat' u 1:3 w l
replot 'drast070.dat' u 1:4 w l
unset multiplot

set term x11 persist 3
set multiplot
unset autoscale
set xrange [0:400]
set yrange [0:3.2]
plot 'drast090.dat' u 1:2 w l
replot 'drast090.dat' u 1:3 w l
replot 'drast090.dat' u 1:4 w l
unset multiplot

set term x11 persist 3
set multiplot
unset autoscale
set xrange [0:400]
set yrange [0:3.2]
plot 'drast200.dat' u 1:2 w l
replot 'drast200.dat' u 1:3 w l
replot 'drast200.dat' u 1:4 w l
unset multiplot