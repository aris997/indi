splot '< tail -n 40000 rho100.dat' u 2:3:4 w l // con questo si plotta la coda del file

//plotta dalla riga 2000 alla riga 4000
plot 'costante.dat' every ::2000::4000 w l
