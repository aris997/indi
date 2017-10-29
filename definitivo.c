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

double ask(char* argomento, double, double);


int main(int argc, char *argv[]) {

	int k, kstop, l;
	long int i, steps, Hstep;
	double x0, y0, z0;
  double C, Co, tmax, delta_periodo;
  //char filename[30];
	
  double *dperiodo;

	FILE *output;
	//FILE *input;
	FILE *periodo;
	FILE *cost;
  //FILE *costerr;
	//FILE *piripicchio;
  //FILE *rho_variable;
  //FILE *drastico;
  //FILE *gpscript;
	
  vector p, pold;
  vector S, Sold;
  vector e;
	pars c;

    c.a = ask("a", 0, 10);
    c.b = ask("b", 0, 10);
    c.rho = ask("rho", 0, 5);
    p.x = ask("x0", 0, 10);
    p.y = ask("y0", 0, 10);
    p.z = ask("z0", 0, 10);
    c.dt = ask("dt", 0.0001, 1.);
    tmax = ask("tmax", 50., 20000.);

	
	output = fopen("output.dat", "w");
	fprintf(output, "#t\t\t\tx\t\t\t\t\ty\t\t\t\t\tz\n");

	//periodo = fopen("periodo.dat", "w");
	//fprintf(periodo, "#periodo\n");

	//costerr = fopen("costante_dt.dat", "w");
  cost = fopen("cost.dat", "w");
  //periodo = fopen("periodo_pp1.dat", "w");
  //fprintf(periodo, "#periodo\n");

	x0 = p.x; y0 = p.y;	z0 = p.z;
  Co = log(x0) - x0 + log(y0) - y0;

  k=0; kstop = NUM_PERIODI; j=0;

	steps = (long int)(tmax/c.dt);
  Hstep = steps/2;

  dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
  if ( dperiodo == NULL ) {printf("malloc error\n"); exit(-1);}
  else dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));


	for (i=0; i<=steps; i++) {

		pold = strcopy(p);
		p = RK4(p, c);
		
    //fprintf(output, "%.8lf\t%.14lf\t%.14lf\t%.14lf\n", ((double)i+1)*c.dt, p.x, p.y, p.z);
		C = log(p.x) - p.x + log(p.y) - p.y + p.z;
		//fprintf(cost, "%lf %.1lf %.18lf\n", ((double)i+1)*c.dt, C, C/Co - 1);
    //printf("%.14lf %.14lf\n", C, Co);

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
          delta_periodo = dperiodo[k-1] - dperiodo[k-2];
          //if (k > 1) fprintf(periodo, "%.8lf %.10e\n", delta_periodo, delta_periodo/c.dt);
          
          if (k == (kstop - 4) ) { 
            dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); 
            kstop *= 2; printf("#reallocated dperiod memory\n");
          }
        }
  		}
    }
  }

  free(dperiodo);
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
    //fprintf(costerr, "%.14lf %.24lf\n", c.dt, fabs(C/Co - 1) );
  }

      /****************************/
      /****************************/
      /****************************/
      /****INIZIO SECONDA PARTE****/printf("#Parte 2\n");
      /****************************/
      /****************************/
      /****************************/

  printf("#Se rho massimo e minimo sono uguali si eseguira un solo ciclo\n");



  double Mrho = ask("rho massimo", 0, 5);
  double mrho = ask("rho minimo", 0, Mrho);
  double drho, passo;
  double dperiodototale;
  int w;

  if (Mrho != mrho ) {
    passo = ask("passo, non troppo piccolo per intervalli grandi", 0.001, fabs(Mrho-mrho));
    drho = fabs(Mrho-mrho)*passo;
  }
  else {
    drho = 1.;
  }



	for (c.rho=mrho; c.rho<=Mrho; c.rho+=drho) {


    dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double));
    if (dperiodo == NULL ) { printf("malloc error\n"); exit(-1); }
    dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));
    kstop = NUM_PERIODI;


    c.a = 0.5;
    c.b = 0.1;
    x0 = 1.2;
    y0 = 0.5;
    z0 = 0.3;
    c.dt = 0.01;
    tmax = 1000.;


		//sprintf(filename, "rho%0*d.dat", width, (int)(c.rho*100.));
		//rho_variable = fopen(filename, "w");

		p.x = x0; p.y = y0; p.z = z0;
    k = 0;

    steps = (long int)(tmax/c.dt);
    Hstep = steps/2;

		for (i=0; i <= steps; i++) {

      pold = strcopy(p);
      p = RK4(p, c);
      
      if ( i == Hstep ) { 
        S = strcopy(p); 
        Sold = strcopy(pold); 
        e = intorno(S, Sold, 1);
      }
			
      else if ( i > Hstep/2. ) { 

        if (e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {
          //controllo se la direzione del passaggio è la stessa
          if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
            
            dperiodo[k] = (double)i*c.dt;
            k++; //aggiornamneto contatori

            if (k == (kstop - 4) ) { //
              dperiodo = (double *) realloc( dperiodo, 2 * kstop * sizeof(double)); 
              kstop *= 2; 
              //printf("#reallocated dperiod memory\n"); 
            }
          }
        }
      }
    }
    if (k > 0) {
      for (w=0; w<(k-1); w++) dperiodototale += dperiodo[w+1]-dperiodo[w];
      dperiodototale/=(k-1);
      printf("#periodo rho: %1.4lf \tT %.8lf \tincontri %d\n", c.rho, dperiodototale, k+1);
    }
    else {
      printf("#periodo NON esiste per rho: %lf\n", c.rho);
    }


    free(dperiodo);
    //fclose(rho_variable);
    //fclose(periodo);

	}


  //fclose(costerr);
	
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

struct vector intorno(vector S, vector Sold, double i) {
  vector e;

  e.x = fabs(i*(S.x-Sold.x)/2.);
  e.y = fabs(i*(S.y-Sold.y)/2.);
  e.z = fabs(i*(S.z-Sold.z)/2.);

  return e;
}

double ask(char* argomento, double a, double b) {
  double w;

  printf("#Inserire %s\n", argomento);

  do {
    scanf("%lf", &w);
    if ( w < a || w > b ) {
      printf("#Valore errato, %s assumere valori [%.2lf:%.2lf]\n", argomento, a, b);
    }
  } while ( w < a || w > b );

  return w;
}

