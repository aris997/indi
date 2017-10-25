#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NUM_PERIODI 256

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

double f(double, double, double, pars);
double g(double, double, double, pars);
double h(double, double, double, pars);



int main() {

  int tmax, k=0, j=0;
  long int i, steps;
  double x0, y0, z0, C;
  double dperiodo[256];

  FILE *output;
  FILE *input;
  FILE *periodo;
  FILE *cost;

  vector p;
  pars c;

  input = fopen("input.dat", "r");
  fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &p.x, &p.y, &p.z, &c.dt, &tmax);

  output = fopen("output.dat", "w");
  fprintf(output, "#t x y z\n");

  periodo = fopen("perio.dat", "w");
  fprintf(periodo, "#periodo ");

  cost = fopen("costante.dat", "w");
  fprintf(cost, "#Costante\n");

  x0 = p.x;
  y0 = p.y;
  z0 = p.z;

  fprintf(output, "%.8lf %.8lf %.8lf %.8lf\n", 0., x0, y0, z0);


  steps = (tmax/c.dt);

	for (i=0; i<steps; i++) {


    p = RK4(p, c);
    fprintf(output, "%.8lf %.16lf %.16lf %.16lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);

    C = log(p.x) - p.x + log(p.y) - p.y + p.z;
    fprintf(cost, "%.14lf\n", C);

    if ( (double)i*dt == 10.) {
      
    }
	}
  

  fprintf(periodo, "%d\n", k-1);
  for (i=0; i<k-1; i++) {
    fprintf(periodo, "%lf\n", dperiodo[i+1]-dperiodo[i]);
  }





  FILE *rho_variable;
  char filename[30];

  steps = (long int)(tmax/c.dt);


  for (j=0; j<2; j++) {

    fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &x0, &y0, &z0, &c.dt, &tmax);
    sprintf(filename, "rho%d.dat", (int)(c.rho*100.));

    printf("%d ", (int)(c.rho*100.));

    rho_variable = fopen(filename, "w");

    p.x = x0;
    p.y = y0;
    p.z = z0;

    for (i=0; i<steps; i++) {
      p = RK4(p, c);
      fprintf(rho_variable, "%.4lf %.8lf %.8lf %.8lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);
    }

  fclose(rho_variable);

  }


  fclose(periodo);
  fclose(cost);
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