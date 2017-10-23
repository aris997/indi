#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct vector { 	//definisco un typedef vector per un vettore nello spazio.
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

struct vector RK4(vector n, pars c, double t);

double f(double z, double t, vector n, pars c);
double g(double z, double t, vector n, pars c);
double h(double z, double t, vector n, pars c);

int main() {
	int tmax; 	          // definisco i parametri
	long int i, steps;		// scelgo un long int per la variabile del ciclo
	FILE *output;			    // def un puntatore per il file out

  //per semplificare abbrevio i nomi delle struct	
  vector p; //p sta per posizione
  pars c;  //c per costanti

  double x0 = 1.2;
  double y0 = 0.5;
  double z0 = 0.;

  c.a = 0;
  c.b = 0;
  c.rho = 0.2;
  c.dt = 0.001;

  tmax = 100;


  output = fopen("output.dat", "w"); //iniz il file
  fprintf(output, "#t x y z\n");
  fprintf(output, "%.8lf %.8lf %.8lf %.8lf\n", 0., x0, y0, z0);

  //il vector p assume le variabili iniziali
  p.x = x0;
  p.y = y0;
  p.z = z0;
  

	steps = (tmax/c.dt); //assegno a steps il numero di passi di integrazione

	for (i=0; i<steps; i++) {  
   
    p = RK4(p, c, (double)(i + 1) * c.dt);
    fprintf(output, "%.4lf %.8lf %.8lf %.8lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);

	}

  fclose(output);
	exit(0);
}

struct vector RK4(vector n, pars c, double t) { //chiamo n il vettore così n mi richiama alla mente a che passo di int sono
  //i passi xi, yi e zi sono le velocità
  
  double x1, x2, x3, x4;  //definisco i passi di RK4 temporanee per x
  double dt = c.dt;

  x1 = f(n.x, t, n, c) * dt;
  x2 = f(n.x + 0.5 * x1, t + 0.5 * dt, n, c) * dt;
  x3 = f(n.x + 0.5 * x2, t + 0.5 * dt, n, c) * dt;
  x4 = f(n.x + x3, t + dt, n, c) * dt;

  n.x += (x1 + 2.*x2 + 2.*x3 + x4)/6.;

  double y1, y2, y3, y4;  //definisco i passi di RK4 temporanee per y

  y1 = g(n.y, t, n, c) * dt;
  y2 = g(n.y + 0.5 * y1, t + 0.5 * dt, n, c) * dt;
  y3 = g(n.y + 0.5 * y2, t + 0.5 * dt, n, c) * dt;
  y4 = g(n.y + y3, t + dt, n, c) * dt;

  n.y += (y1 + 2.*y2 + 2.*y3 + y4)/6.;

  double z1, z2, z3, z4;  //definisco i passi di RK4 temporanee per z

  z1 = h(n.z, t, n, c) * dt;
  z2 = h(n.z + 0.5 * z1, t + 0.5 * dt, n, c) * dt;
  z3 = h(n.z + 0.5 * z2, t + 0.5 * dt, n, c) * dt;
  z4 = h(n.z + z3, t + dt, n, c) * dt;

  n.y += (z1 + 2.*z2 + 2.*z3 + z4)/6.;

  return n;
}


/****LE FUNZIONI f, g e h del problema****/

double f(double x, double t, vector n, pars c) {
  return x * ( c.a * (1 - n.x) + 0.5 * (1 - n.y) + c.b * (1 - n.z) );
}

double g(double y, double t, vector n, pars c) {
  return y * ( -0.5 * (1 - n.x) + c.b * (n.y - n.z) );
}

double h(double z, double t, vector n, pars c) {
  return z * ( c.rho * (1 - n.x) + 0.1 * (2. - n.y - n.z) );
}