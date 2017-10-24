#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NUM_PERIODI 128

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

struct vector RK4(vector, pars, double);
struct vector strcopy(vector);

double f(double, double, vector, pars);
double g(double, double, vector, pars);
double h(double, double, vector, pars);

int main() {
  int tmax; 	          // definisco i parametri
  int k=0, j=0;         // k è una variabile necessaria per posizionare i tempi in un array malloc, j mi aiuta nell'ordine
  long int i, steps;		// scelgo un long int per la variabile del ciclo
  double x0, y0, z0, dt;		// 
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
  vector pold; //pold mi conserva i dati ad n-1 per poter calcolare il periodo
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

	p = RK4(p, c, (double)(i + 1) * c.dt);
	fprintf(output, "%.8lf %.16lf %.16lf %.16lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);


  /****Zona del controllo periodo****/
	//if (p.x <= pold.x && p.y - y0 <= y0*0.001 && p.z - z0 <= z0*0.001) { //se periodica dovrebbe tornare al punto di partenza
	if (p.x >= x0*(1-dt*dt) && pold.x <= x0*(1+dt*dt) && p.y >= y0*(1-dt*dt) && pold.y <= y0*(1+dt*dt) && p.z >= z0*(1-dt*dt) && pold.z <= z0*(1+dt*dt)) {
	  dperiodo[k] = c.dt*((double)(i+1)); //temporaneamente metto i tempi in un array
	  k++;  //conto quanti ne salvo
	  /****REALLOC ZONE - in fase di costruzione****/
	  // if (k == NUM_PERIODI * (j+1) - 1) { //k+1, prima di finire completamente la memoria ne assegno della nuova
	  //   dperiodo = realloc(dperiodo, 2 * NUM_PERIODI * sizeof(double)); //raddoppio lo spazio riservato a dperiodo
	  //   j++; //j mi aiuta a tenere ordinata la memoria
	  // }
	}
  
  /****SECONDO PUNTO PARTE 1****/
	fprintf(energy, "%.8lf %.16lf\n", c.dt*((double)(i+1)), (log(p.x) - p.x + log(p.y) - p.y + p.z));


	}
  fclose(output);







/****CALCOLO E SALVATAGGIO DEL PERIODO****/
  fprintf(periodo, "%d\n", k-1);
  
  for (i=0; i<k-1; i++) {
	fprintf(periodo, "%lf\n", dperiodo[i+1]-dperiodo[i]);
  }


/****CONCLUSIONI****/

  free(dperiodo); //heap libero
  fclose(periodo);
  fclose(energy);

	exit(0);
}



struct vector RK4(vector n, pars c, double t) { //chiamo n il vettore così n mi richiama alla mente a che passo di int sono
  //i passi xi, yi e zi sono le velocità
  
  double x1, x2, x3, x4;  //definisco i passi di RK4 temporanee per x
  double dt = c.dt;

  vector p1, p2, p3, p4;

  // x1 = f(n.x, n, c) * dt;
  // y1 = g(n.y, n, c) * dt;
  // z1 = h(n.z, n, c) * dt;

  p1 = q(n.x, n.y, n.z, c) * dt;

  // x2 = f(n.x + 0.5 * x1, n, c) * dt;    //repetita iuvant: i passi vengono fatti per tutte e tre le dimensioni
  // y2 = g(n.y + 0.5 * y1, n, c) * dt;
  // z2 = h(n.z + 0.5 * z1, n, c) * dt;

  p2 = q(n.x + p1.x/2., n.y + p1.y/2., n.z + p1.z/2., c) * dt;

  double y1, y2, y3, y4;  //definisco i passi di RK4 temporanee per y

  // x3 = f(n.x + 0.5 * x2, n, c) * dt;
  // y3 = g(n.y + 0.5 * y2, n, c) * dt;
  // z3 = h(n.z + 0.5 * z2, n, c) * dt;

  p3 = q(n.x + p2.x/2., n.y + p2.y/2., n.z + p2.z/2., c) * dt;

  double z1, z2, z3, z4;  //definisco i passi di RK4 temporanee per z

  // x4 = f(n.x + x3, n, c) * dt;
  // y4 = g(n.y + y3, n, c) * dt;
  // z4 = h(n.z + z3, n, c) * dt;

  p4 = q(n.x + p3.x, n.y + p3.y, n.z + p3.z, c) * dt;

  n.x += (p1.x + 2.*p2.x + 2.*p3.x + p4.x)/6.;
  n.y += (p1.y + 2.*p2.y + 2.*p3.y + p4.y)/6.;
  n.z += (p1.z + 2.*p2.z + 2.*p3.z + p4.z)/6.;

  return n;
}



/****LE FUNZIONI f, g e h del problema(in ordine x, y, z)****/

double f(double x, vector n, pars c) {
  return x * ( c.a * (1 - n.x) + 0.5 * (1 - n.y) + c.b * (1 - n.z) );
}

double g(double y, vector n, pars c) {
  return y * ( -0.5 * (1 - n.x) + c.b * (n.y - n.z) );
}

double h(double z, vector n, pars c) {
  return z * ( c.rho * (1 - n.x) + 0.1 * (2. - n.y - n.z) );
}

struct vector q(double x, double y, double z, pars c) {
  vector p;

  p.x = x * ( c.a * (1 - x) + 0.5 * (1 - y) + c.b * (1 - z) );
  p.y = y * ( -0.5 * (1 - x) + c.b * (y - z) );
  p.z = z * ( c.rho * (1 - x) + 0.1 * (2. - y - z) );

  return p;
}

/****Utili****/
struct vector strcopy(vector a) {
  
  vector new;
  
  new.x = a.x;
  new.y = a.y;
  new.z = a.z;

  return new;
}