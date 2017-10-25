set term x11 persist //questo sceglie il terminale (nel caso di mac aqua o xorg quello che usate)

set xrange [9.95:10.05] //con questo date un range (basta fare yrange o zrange)

set multiplot	//attiva la modalità per sovrapporre grafici (se il range è definito non ci sono problemi di sovrapposizione)


plot 'output.dat' u 1:2 w p //plot regolare con utilizzo delle colonne 1 e 2 with lines

//in alternativa si possono scegliere i punti con il lessico "w p" (with points)

unset multiplot //disattiva la modalità multiplot

replot //si può usare sia con splot che con plot

plot '< tail -n 20000 out.dat' u 2:3  w l //prende le ultime 20000 righe del file out.dat e plotta la seconda e terza colonna con delle linee (non sono sicuro funzioni su mac, vale la pena tentare)

plot 'output.dat' u 1:2 with points