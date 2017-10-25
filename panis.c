#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NUM_PERIODI 256

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

struct vector RK4(vector, pars);
struct vector strcopy(vector);

double f(double, double, double, pars);
double g(double, double, double, pars);
double h(double, double, double, pars);

int main() {
  int tmax; 	          // definisco i parametri
  int k=0;              // k è una variabile necessaria per posizionare i tempi in un array malloc, j mi aiuta nell'ordine
  long int i, steps;		// scelgo un long int per la variabile del ciclo
  double x0, y0, z0, dt, C;		// 
  FILE *output;			    // def un puntatore per il file out
  FILE *input;
  FILE *periodo;
  FILE *energy;

  double *dperiodo; //definisco delta dei periodi

  dperiodo = (double *) malloc(NUM_PERIODI * sizeof(double)); //inizializzo
  if (dperiodo == NULL) { //controllo
	printf("malloc dperiodo failed\n");
	exit(EXIT_FAILURE);
  }
  else {  //azzero
	  dperiodo = (double *) calloc(NUM_PERIODI, sizeof(double));
  }

  //per semplificare abbrevio i nomi delle struct	
  vector p; //p sta per posizione
  vector pold;
  pars c;  //c sta per costanti


/****Raccolgo da input.dat le condizioni iniziali****/
  input = fopen("input.dat", "r");
  fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &p.x, &p.y, &p.z, &c.dt, &tmax);
  fclose(input);

/****Creo un output.dat e scrivo l'ordine dei dati****/
  output = fopen("output.dat", "w"); //iniz il file
  fprintf(output, "#t x y z\n");

/****Apro un altro file per i periodi****/
  periodo = fopen("perio.dat", "w");
  fprintf(periodo, "#periodo ");

/**Apro un file per salvare l'energy****/
  energy = fopen("energy.dat", "w");
  fprintf(energy, "#C\n");


  //x0, y0 e z0
  x0 = p.x;
  y0 = p.y;
  z0 = p.z;

  fprintf(output, "%.8lf %.8lf %.8lf %.8lf\n", 0., x0, y0, z0);


  steps = (tmax/c.dt); //assegno a steps il numero di passi di integrazione
  dt = c.dt;

/****Zona principale del codice****/
	for (i=0; i<steps; i++) {
    pold = strcopy(p);  //ricopio il vettore dentro uno nuovo per poter controllare il periodo
    
    p = RK4(p, c);
    fprintf(output, "%.8lf %.16lf %.16lf %.16lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);


  /****Zona del controllo periodo****/
  //if ()
  
  /****SECONDO PUNTO PARTE 1****/
    C = log(p.x) - p.x + log(p.y) - p.y + p.z;
    fprintf(energy, "%.24lf\n", C);


	}
  
  fclose(output);


  FILE *cosettero;

  cosettero = fopen("cosettero.dat", "w");
  fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &p.x, &p.y, &p.z, &c.dt, &tmax);

  steps = (long int)(tmax/c.dt);

  for (c.rho=0.75; c.rho<=1.00; c.rho+=0.05) {

    for (i=0; i<steps; i++) {
      p = RK4(p, c);
      fprintf(cosettero, "%.4lf %.8lf %.8lf %.8lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);
    }
  }






/****CALCOLO E SALVATAGGIO DEL PERIODO****/
  fprintf(periodo, "%d\n", k-1);
  
  for (i=0; i<k-1; i++) {
	fprintf(periodo, "%lf\n", dperiodo[i+1]-dperiodo[i]);
  }


/****CONCLUSIONI****/
  fclose(cosettero);
  free(dperiodo); //heap libero
  fclose(periodo);
  fclose(energy);

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



/****LE FUNZIONI f, g e h del problema(in ordine x, y, z)****/

double f(double x, double y, double z, pars c) {
  return x * ( c.a * (1 - x) + 0.5 * (1 - y) + c.b * (1 - z) );
}

double g(double x, double y, double z, pars c) {
  return y * ( -0.5 * (1 - x) + c.b * (y - z) );
}

double h(double x, double y, double z, pars c) {
  return z * ( c.rho * (1 - x) + 0.1 * (2. - y - z) );
}




/****Utili****/
struct vector strcopy(vector a) {
  
  vector new;
  
  new.x = a.x;
  new.y = a.y;
  new.z = a.z;

  return new;
}