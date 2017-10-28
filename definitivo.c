#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NUM_PERIODI 16

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
struct vector intorno(vector, vector, double);

double f(double, double, double, pars);
double g(double, double, double, pars);
double h(double, double, double, pars);


int main(int argc, char *argv[]) {

	int k, kstop, j, width, l;
	long int i, steps, Hstep;
	double x0, y0, z0;
  double C, Co, tmax;
  //char filename[30];
	
  double *dperiodo;

	FILE *output;
	FILE *input;
	FILE *periodo;
	FILE *cost;
  FILE *costerr;
	//FILE *piripicchio;
  //FILE *rho_variable;
  //FILE *drastico;
  //FILE *gpscript;
	
  vector p, pold;
  vector S, Sold;
  vector e;
	pars c;

  if ( argc == 10 ) {
    c.a = atof(argv[2]);
    c.b = atof(argv[3]);
    c.rho = atof(argv[4]);
    p.x = atof(argv[5]);
    p.y = atof(argv[6]);
    p.z = atof(argv[7]);
    c.dt = atof(argv[8]);
    tmax = atof(argv[9]);
  
    if (tmax < 50. && c.dt >= 1.) {
      printf("#tmax troppo piccolo o dt troppo grande\nEsecuzione invalidata\n");
      exit(-1);
    }
  }



  else {
    printf("#Inserire (-1) per presa dati da file input\n");
    printf("#Inserire: a b rho x0 y0 z0 dt tmax da linea di comando per usare tali valori scelti\n");
    printf("#l esecuzione prenderà i dati di default\n");
    p.x = 1.2;
    p.y = 0.5;
    p.z = 0.;
    c.a = 0;
    c.b = 0;
    c.rho = 0.2;
    c.dt = 0.01;
    tmax = 50.;
  }
	
	output = fopen("output.dat", "w");
	fprintf(output, "#t\t\t\tx\t\t\t\t\ty\t\t\t\t\tz\n");

	//periodo = fopen("periodo.dat", "w");
	//fprintf(periodo, "#periodo\n");

	costerr = fopen("costante_dt.dat", "w");
  cost = fopen("cost.dat", "w");


	x0 = p.x; y0 = p.y;	z0 = p.z;
  Co = log(x0) - x0 + log(y0) - y0;

  k=0; kstop = NUM_PERIODI; j=0;

	steps = (long int)(tmax/c.dt);
  Hstep = steps/2;

  //dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
  //if ( dperiodo == NULL ) {printf("malloc error\n"); exit(-1);}
  //else dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));


	for (i=0; i<=steps; i++) {

		pold = strcopy(p);
		p = RK4(p, c);
		fprintf(output, "%.8lf\t%.14lf\t%.14lf\t%.14lf\n", ((double)i+1)*c.dt, p.x, p.y, p.z);
		C = log(p.x) - p.x + log(p.y) - p.y + p.z;
		fprintf(cost, "%lf %.14lf\n", ((double)i+1)*c.dt, C/Co -1);
    //printf("%.14lf %.14lf\n", C, Co);
/*
    if ( i == Hstep ) { 
      S = strcopy(p); 
      Sold = strcopy(pold); 
      e = intorno(S, Sold, 1);
    }

    if ( i > Hstep ) {
  		if ( e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {
  			
        if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
          
          dperiodo[k] = (double)i*c.dt;
          k++;
          printf("#%d\n", k);
          
          if (k > 1) fprintf(periodo, "%lf %.8lf\n", (double)i*c.dt, dperiodo[k-1] - dperiodo[k-2]);
          
          if (k == (kstop - 4) ) { 
            dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); 
            kstop *= 2; printf("#reallocated dperiod memory\n");
          }
        }
  		}
    }*/
  }

  //free(dperiodo);
  fclose(cost);
  fclose(output);


  printf("#Vario dt per osservare l'andamento di C/Co - 1\n");
  
  /****ricalcolo il moto fino a 10.t variando il dt****/
  for (c.dt=1.; c.dt>=0.0001; c.dt/=10.) {
    tmax = 10.;
    steps = (long int)(tmax/c.dt);
    p.x = x0;
    p.y = y0;
    p.z = z0;
    //printf("%lf\n", c.dt);
    for (i=0; i<=steps; i++) { p = RK4(p, c); }
    C = log(p.x) - p.x + log(p.y) - p.y + p.z;
    //printf("%.14lf %.14lf %.14lf %.24lf\n", p.x, p.y, p.z, C);
    fprintf(costerr, "%.14lf %.24lf\n", c.dt, fabs(C/Co - 1) );
  }

      /****************************/
      /****************************/
      /****************************/
      /****INIZIO SECONDA PARTE****/printf("#INIZIO SECONDA PARTE\n");
      /****************************/
      /****************************/
      /****************************/
/*
	for (j=0; j<2; j++) {

    //printf("#E %d\n", j+1);

    dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
    if (dperiodo == NULL ) { printf("malloc error\n"); exit(-1); }
    dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));

		//fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &c.a, &c.b, &c.rho, &x0, &y0, &z0, &c.dt, &tmax);

		sprintf(filename, "rho%0*d.dat", width, (int)(c.rho*100.));
		rho_variable = fopen(filename, "w");

		p.x = x0; p.y = y0; p.z = z0;
    k = 0; l = 0;

    steps = (long int)(tmax/c.dt);
    Hstep = steps/2;

    sprintf(filename, "periodorho%0*d.dat", width, (int)(c.rho*100.));
    periodo = fopen(filename, "w");
    fprintf(periodo, "#Condizioni Iniziali\n");
    fprintf(periodo, "#%lf %lf %lf %lf %lf %lf %lf %lf\n", p.x, p.y, p.z, c.a, c.b, c.rho, c.dt, tmax);

		for (i=0; i <= steps; i++) {

      pold = strcopy(p);
      p = RK4(p, c);
      
      if ( i == Hstep ) { 
        S = strcopy(p); 
        Sold = strcopy(pold); 
        e = intorno(S, Sold, 1);
      }
			
      else if ( i > Hstep/2. ) { 
        fprintf(rho_variable, "%.4lf %.8lf %.8lf %.8lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);

        // controllo se il moto del punto si trova nei punti scelti S e Sold 
        if (e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {

          //controllo se la direzione del passaggio è la stessa
          if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
            
            dperiodo[k] = (double)i*c.dt;
            k++; l++; //aggiornamneto contatori

            if (k == (kstop - 4) ) { //
              dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); 
              kstop *= 2; 
              printf("#reallocated dperiod memory\n"); 
            }
            
            if (k > 1) { 
              fprintf(periodo, "%.8lf %.14lf\n", (double)i*c.dt, dperiodo[k-1]-dperiodo[k-2]); 
            }
            
            // if ( k > 2 && dperiodo[k-1]-dperiodo[k-2] < (dperiodo[k-2]-dperiodo[k-3])/4.) {
            //   printf("#A S I N T O T O   T I M E\t %lf\n",( ( (double)i*c.dt) - 3 ) );
            //   fprintf(rho_variable, "#Asintotico a {x:%.8lf, y:%.8lf, z:%.8lf, t:%lf}\n", p.x, p.y, p.z, (double)i*c.dt);
            //   i = steps;
            // }
          }
        }
      }
    }



    free(dperiodo);
    fclose(rho_variable);
    fclose(periodo);

	}

*/
  printf("sasd\n");

/****RIORDINARE QUESTA SEZIONE (RIGUARDA LA PARTE 2)****/
/*
	c.a = 0.5;
	c.b = 0.1;
	x0 = 1.2;
	y0 = 0.5;
	z0 = 0.3;
	c.dt = 0.01;
	tmax = 200;

	steps = (long int)(tmax/c.dt);
  Hstep = steps/2;


	for (j=0; j<=50; j++)  { //100 passi per arrivare da 1.2 a 1.5 rho


		c.rho = 1.200;
		c.rho += 0.006*(double)j; //incremento dell'1% dell'intervallo
		p.x = x0;
		p.y = y0;
		p.z = z0;
    k = 0;

    
    dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
    if (dperiodo == NULL ) { printf("malloc error\n"); exit(-1); }
    dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));
    kstop = NUM_PERIODI;
		

    for (i = 0; i<=steps; i++) { //steps dovrebbe essere 8000

			pold = strcopy(p);
			p = RK4(p, c);

			if (i == Hstep) {
        e = intorno(S, Sold, 20);
        if (e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {
          
          if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
            //printf("#cambiamento non drastico con rho:%.8lf\n", c.rho);
          }
        }
        else {
          //printf("#cambiamento drastico, con rho:%.8lf\n", c.rho);
        }      


        Sold = strcopy(pold);
        S = strcopy(p);
        e = intorno(S, Sold, 0.7);
      }
			
      if (i > Hstep) {

        //printf("%lf %lf %lf ", p.x, p.y, p.z);
				
        if (e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {
          if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
            //printf("#periodo %lf\n", (double)i*c.dt); 
            dperiodo[k] = (double)i*c.dt; k++;
        
            //if (k > 1) printf("%lf %.8lf ", c.rho, dperiodo[k-1]-dperiodo[k-2]);
           
            if (k == (kstop - 4) ) { //
              dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); 
              kstop *= 2; 
              //printf("#reallocated dperiod memory\n"); 
            }
          }
				}
        //printf("\n");
			}
		}
  free(dperiodo);
	}
*/
  printf("porcodio");
  
  fclose(costerr);
	//fclose(input);

	exit(0);
}//////////////////////////////// end main



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

struct vector intorno(vector S, vector Sold, double i) {
  vector e;

  e.x = fabs(i*(S.x-Sold.x)/2.);
  e.y = fabs(i*(S.y-Sold.y)/2.);
  e.z = fabs(i*(S.z-Sold.z)/2.);

  return e;
}