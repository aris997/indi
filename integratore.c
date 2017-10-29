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

int main () {
  
  long int passaggi, stoppassaggi;
  long int i, k, w;
  long int passi, metapassi;
  long int MAX;

  //double C, Co;
  double x0, y0, z0;
  double tmax;
  double periodosomma;

  double *periodo;

  vector p, pold;
  vector S, Sold;
  vector e;
  pars c;
  
  FILE *periodifile;
  FILE *output;
  //FILE *errore;

  c.a = ask("a", 0, 10);
  c.b = ask("b", 0, 10);
  p.x = ask("x0", 0, 10);
  p.y = ask("y0", 0, 10);
  p.z = ask("z0", 0, 10);
  c.dt = ask("dt", 0.0001, 1.);
  tmax = ask("tmax", 50., 20000.);

  printf("#Se rho massimo e minimo sono uguali si eseguira un solo ciclo\n");

  double Mrho = ask("rho massimo", 0, 5);
  double mrho = ask("rho minimo", 0, Mrho);
  double rhopasso, passo;

  if (Mrho != mrho ) {
    passo = ask("passo, non troppo piccolo per intervalli grandi", 0.001, fabs(Mrho-mrho));
    rhopasso = fabs(Mrho-mrho)*passo;
    MAX = (long int)(1/passo);
  }
  else {
    MAX = 0;
  }


  x0 = p.x; y0 = p.y; z0 = p.z;
  //Co = log(x0) - x0 + log(y0) - y0;

  periodifile = fopen("periodi.dat", "w");
  output = fopen("output.dat", "w");

  for (w=0; w<=MAX; w++) {

    c.rho = mrho + ((double)(w))*rhopasso;
    //printf("#ciclo numero %ld \n", w+1);
    
    periodo = (double *) malloc(NUM_PERIODI * sizeof(double));
    if (periodo == NULL ) { printf("malloc error\n"); exit(-1); }
    periodo = (double *) calloc(NUM_PERIODI, sizeof(double));
    stoppassaggi = NUM_PERIODI;

    p.x = x0; p.y = y0; p.z = z0;
    passaggi = 0;

    passi = (long int)(tmax/c.dt);
    metapassi = passi/2;

    for (i=0; i <= passi; i++) {

      pold = strcopy(p);
      p = RK4(p, c);
      //fprintf(output, "%lf %.14lf %.14lf %.14lf\n", (double)i*c.dt, p.x, p.y, p.z);
      
      if ( i == metapassi ) { 
        S = strcopy(p); 
        Sold = strcopy(pold); 
        e = intorno(S, Sold, 1);
      }
      
      else if ( i > metapassi ) { 

        if (e.x >= fabs(p.x - S.x) && e.y >= fabs(p.y - S.y) && e.z >= fabs(p.z - S.z) ) {
          //controllo se la direzione del passaggio è la stessa
          if (e.x >= fabs(pold.x - Sold.x) && e.y >= fabs(pold.y - Sold.y) && e.z >= fabs(pold.z - Sold.z) ) {
            
            periodo[passaggi] = (double)i*c.dt;
            passaggi++; //aggiornamneto contatori

            if (passaggi == (stoppassaggi - 4) ) { //
              periodo = (double *) realloc( periodo, 2 * stoppassaggi * sizeof(double)); 
              stoppassaggi *= 2; 
              //printf("#reallocated periodo memory\n"); 
            }
          }
        }
      }
    }

    for (k=0; k<(passaggi-1); k++) periodosomma += periodo[k+1]-periodo[k];
    periodosomma/=(double)(passaggi-1);

    if (passaggi > 0 && periodosomma > c.dt*2) {
      printf("#periodo rho: %1.4lf \tT %.8lf \tincontri %ld\n", c.rho, periodosomma, passaggi+1);
      fprintf(periodifile,"%.4lf \t%.8lf \t %ld\n", c.rho, periodosomma, passaggi+1);
    }
    
    else {
      printf("#Asintoto esiste per rho: %lf\n", c.rho);
      fprintf(periodifile, "%.4lf \t0.0 \t0\n", c.rho);
    }


    free(periodo);

  }



  //errore = fopen("errore.dat", "w");
/*
  for (c.dt=1.; c.dt>=0.0001; c.dt/=10.) {
    tmax = 10.;
    passi = (long int)(tmax/c.dt);
    p.x = x0;
    p.y = y0;
    p.z = z0;
    //printf("%lf\n", c.dt);
    for (i=0; i<=passi; i++) { p = RK4(p, c); }
    C = log(p.x) - p.x + log(p.y) - p.y + p.z;
    //printf("%.14lf %.14lf %.14lf %.24lf\n", p.x, p.y, p.z, C);
    //fprintf(errore, "%.14lf %.24lf\n", c.dt, fabs(C/Co - 1) );
  }
*/
  fclose(periodifile);
  fclose(output);
  //fclose(errore);

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