#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main () {

	double *periodo;
	periodo = (double *) malloc(100 * sizeof(double));
	if (periodo == NULL ) {printf("malloc error\n"); exit(-1);}
	else periodo = (double *) calloc(100, sizeof(double));
	for (int i=0; i < 101; i++) {printf("%.14lf %d\n", periodo[i], i);}
	periodo = (double *) realloc(periodo, 0.5 * 100 * sizeof(double));
	for (int i=0; i < 51; i++) {printf("%.14lf %d\n", periodo[i], i);}

	return 0;
}