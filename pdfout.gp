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

set output ""