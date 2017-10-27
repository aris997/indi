#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NUM_PERIODI 32

typedef struct vector {
		double x;
		double y;
		double z;
}vector;

typedef struct pars {
	double a;
	double b;
	double rho;
	double dt;
}pars;

struct vector RK4(vector, pars);
struct vector strcopy(vector);
struct vector intorno(vector, vector);

double f(double, double, double, pars);
double g(double, double, double, pars);
double h(double, double, double, pars);


int main(int argv, char *argc[]) {

	int tmax, k, kstop, j, width;
	long int i, steps;
	double x0, y0, z0;
  double C, Co, Hstep;
  char filename[30];
	
  double *dperiodo;

	FILE *output;
	FILE *input;
	FILE *periodo;
	FILE *cost;
  FILE *costerr;
	//FILE *piripicchio;
  FILE *rho_variable;
  //FILE *drastico;
  //FILE *gpscript;
	
  vector p, pold;
  vector S, Sold;
  vector e;
	pars c;

	input = fopen("input.txt", "r");
	fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &p.x, &p.y, &p.z, &c.dt, &tmax);

	output = fopen("output.dat", "w");
	fprintf(output, "#t x y z\n");

	periodo = fopen("periodo.dat", "w");
	//fprintf(periodo, "#periodo\n");

	cost = fopen("costante.dat", "w");
	costerr = fopen("costanterrore.dat", "w");

	//fprintf(cost, "#Costante\n");

	//piripicchio = fopen("piripicchio.dat", "w");

	x0 = p.x; y0 = p.y;	z0 = p.z;

  k=0; kstop = NUM_PERIODI; j=0;

  width = 3;

	fprintf(output, "%.8lf %.8lf %.8lf %.8lf\n", 0., x0, y0, z0);

	Co = log(x0) - x0 + log(y0) - y0;

	steps = (long int)(tmax/c.dt);

  dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
  if ( dperiodo == NULL ) {printf("malloc error\n"); exit(-1);}
  else dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));

  Hstep = steps/2;




	for (i=0; i<steps; i++) {

		pold = strcopy(p);

		p = RK4(p, c);
		//fprintf(output, "%.8lf %.14lf %.14lf %.14lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);

		C = log(p.x) - p.x + log(p.y) - p.y + p.z;
		//fprintf(cost, "%.14lf\n", C);
    //printf("%.14lf %.14lf\n", C, Co);

    if ( i == Hstep ) { S = strcopy(p); Sold = strcopy(pold); e = intorno(S, Sold); }

		if ( e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {
			if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
        dperiodo[k] = (double)i*c.dt;
        k++;
        printf("%d\t", k);
        //printf("%d\t%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", k, Sold.x, S.x, pold.x, p.x, Sold.y, S.y, pold.y, p.y, Sold.z, S.z, pold.z, p.z, (double)i*c.dt);
        if (k > 1) fprintf(periodo, "%lf %.8lf\n", (double)i*c.dt, dperiodo[k-1] - dperiodo[k-2]);
        printf("\n");
//        printf("%d %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf\n", k, Sold.x, S.x, pold.x, p.x, Sold.y, S.y, pold.y, p.y, Sold.z, S.z, pold.z, p.z, (double)i*c.dt);
        if (k == (kstop - 4) ) { dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); kstop *= 2; printf("reallocated dperiod memory\n");}
      }
		}
	}

  free(dperiodo);



/****ricalcolo il moto fino a 10.t variando il dt****/
  for (c.dt=1.; c.dt>=0.0001; c.dt/=2.) {
    tmax = 10;
    steps = (long int)(tmax/c.dt);
    p.x = x0;
    p.y = y0;
    p.z = z0;
    //printf("%lf\n", c.dt);
    for (i=0; i<=steps; i++) { p = RK4(p, c); }
    C = log(p.x) - p.x + log(p.y) - p.y + p.z;
    //printf("%.14lf %.14lf %.14lf %.24lf\n", p.x, p.y, p.z, C);
    fprintf(costerr, "%.14lf %.52lf\n", c.dt, fabs(C/Co - 1) );
  }

      /****************************/
      /****************************/
      /****************************/
      /****INIZIO SECONDA PARTE****/printf("INIZIO SECONDA PARTE\n");
      /****************************/
      /****************************/
      /****************************/

	for (j=0; j<2; j++) {

    dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
    if (dperiodo == NULL ) { printf("malloc error\n"); exit(-1); }
    dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));

		fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &x0, &y0, &z0, &c.dt, &tmax);

		sprintf(filename, "rho%0*d.dat", width, (int)(c.rho*100.));
		rho_variable = fopen(filename, "w");

		p.x = x0; p.y = y0; p.z = z0;
    k = 0; j = 0;

    steps = (long int)(tmax/c.dt);
    Hstep = steps/2;

    sprintf(filename, "periodorho%0*d.dat", width, (int)(c.rho*100.));
    periodo = fopen(filename, "w");
    fprintf(periodo, "#Condizioni Iniziali\n");
    fprintf(periodo, "#%lf %lf %lf %lf %lf %lf %lf %d\n", p.x, p.y, p.z, c.a, c.b, c.rho, c.dt, tmax);

		for (i=0; i <= steps; i++) {

      pold = strcopy(p);
      p = RK4(p, c);
      
      if ( i == Hstep ) { S = strcopy(p); Sold = strcopy(pold); e = intorno(S, Sold); }
			if ( i > Hstep ) { 
        fprintf(rho_variable, "%.4lf %.8lf %.8lf %.8lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);
        
        if (e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {  // controllo se il moto del punto si trova nei punti scelti S e Sold 
          if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
            dperiodo[k] = (double)i*c.dt;
            k++; j++; //aggiornamneto contatori
            if (k == (kstop - 4) ) { dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); kstop *= 2; printf("reallocated dperiod memory\n"); }
            if (k > 1) { fprintf(periodo, "%.8lf %.14lf\n", (double)i*c.dt, dperiodo[k-1]-dperiodo[k-2]); }
          }
        }
        if (j > (3/c.dt) ) {
          printf("A S I N T O T O   T I M E\t %lf\n",( ( (double)i*c.dt) - 3 ) );
          fprintf(rho_variable, "#TROVATO ASINTOTO PER I PUNTI x:%.8lf y:%.8lf z:%.8lf t:%lf\n", p.x, p.y, p.z, (double)i*c.dt);
          i = steps + 1;
        }
        if (i%( (long int)(5/c.dt) ) == 0) j = 0; //ogni 5 secondi controlla quante volte e stato assunto il target
      }
		}

    free(dperiodo);
    fclose(rho_variable);
    fclose(periodo);

	}





	//gpscript = fopen("mscript.gp", "w");
	//fprintf(gpscript, "set term x11 persist\nset multiplot\nunset autoscale\nset yrange [0.5:8]\nset xrange [200:400]\nunset key\nplot");




/****RIORDINARE QUESTA SEZIONE (RIGUARDA LA PARTE 2)****/
  /*
	c.a = 0.5;
	c.b = 0.1;
	x0 = 1.2;
	y0 = 0.5;
	z0 = 0.3;
	c.dt = 0.1;
	tmax = 800;
	k = 0;
	//As = 0;
	e = c.dt/2.;

	steps = (long int)(tmax/c.dt);


	for (j=0; j<=100; j++)  { //100 passi per arrivare da 1.2 a 1.5 rho


		c.rho = 1.197;
		c.rho += 0.003*(double)j; //incremento dell'1% dell'intervallo

		//sprintf(filename, "cambiamento_drast_rho%0*d.dat", width, j);
		// fprintf(gpscript, " '%s' every ::5000::8000 u 1:4 w l\nreplot", filename);
		// drastico = fopen(filename, "w");

		sprintf(filename, "peridi.dat");
		// //completare la scrittura dello script gp
		periodo = fopen(filename, "w");
		fprintf(periodo, "%lf ", c.rho);

		p.x = x0;
		p.y = y0;
		p.z = z0;

		for (i = 0; i<=steps; i++) { //steps dovrebbe essere 8000

			pold = strcopy(p);

			p = RK4(p, c);

			if (i == 4000) { Sold = strcopy(pold); S = strcopy(p); }
			if (i > 4000) {
				//fprintf(drastico, "%.8lf %.14lf %.14lf %.14lf\n", (double)i*c.dt, p.x, p.y, p.z);
				if (e >= fabs(p.x - S.x) && e >= fabs(p.y - S.y) && e >= fabs(p.z - S.z) ) {
          			if (e >= fabs(pold.x - Sold.x) && e >= fabs(pold.y - Sold.y) && e >= fabs(pold.z - Sold.z) ) {
            
          }
				}
			}

			
			//As += asintotico(pold, p, c);
			//if (As == 10000) { i = steps; printf("%lf esiste asintoto\n", c.rho); }
		}

		fclose(drastico);

	}

*/
	fclose(cost);
  fclose(costerr);
	fclose(input);
	fclose(output);

	exit(0);
}



struct vector RK4(vector n, pars c) {
	
	double dt = c.dt;

	vector p1, p2, p3, p4;

	p1.x = f(n.x, n.y, n.z, c) * dt;
	p1.y = g(n.x, n.y, n.z, c) * dt;
	p1.z = h(n.x, n.y, n.z, c) * dt;
	
	p2.x = f(n.x + p1.x/2., n.y + p1.y/2., n.z + p1.z/2., c) * dt;
	p2.y = g(n.x + p1.x/2., n.y + p1.y/2., n.z + p1.z/2., c) * dt;
	p2.z = h(n.x + p1.x/2., n.y + p1.y/2., n.z + p1.z/2., c) * dt;
	
	p3.x = f(n.x + p2.x/2., n.y + p2.y/2., n.z + p2.z/2., c) * dt;
	p3.y = g(n.x + p2.x/2., n.y + p2.y/2., n.z + p2.z/2., c) * dt;
	p3.z = h(n.x + p2.x/2., n.y + p2.y/2., n.z + p2.z/2., c) * dt;
	
	p4.x = f(n.x + p3.x, n.y + p3.y, n.z + p3.z, c) * dt;
	p4.y = g(n.x + p3.x, n.y + p3.y, n.z + p3.z, c) * dt;
	p4.z = h(n.x + p3.x, n.y + p3.y, n.z + p3.z, c) * dt;

	n.x += (p1.x + 2.*p2.x + 2.*p3.x + p4.x)/6.;
	n.y += (p1.y + 2.*p2.y + 2.*p3.y + p4.y)/6.;
	n.z += (p1.z + 2.*p2.z + 2.*p3.z + p4.z)/6.;

	return n;
}



double f(double x, double y, double z, pars c) {
	return x * ( c.a * (1. - x) + 0.5 * (1. - y) + c.b * (1. - z) );
}

double g(double x, double y, double z, pars c) {
	return y * ( -0.5 * (1. - x) + c.b * (y - z) );
}

double h(double x, double y, double z, pars c) {
	return z * ( c.rho * (1. - x) + 0.1 * (2. - y - z) );
}


struct vector strcopy(vector a) {
	
	vector new;
	
	new.x = a.x;
	new.y = a.y;
	new.z = a.z;

	return new;
}

struct vector intorno(vector S, vector Sold) {
  vector e;

  e.x = fabs((S.x-Sold.x)/2.5);
  e.y = fabs((S.y-Sold.y)/2.5);
  e.z = fabs((S.z-Sold.z)/2.5);

  return e;
}