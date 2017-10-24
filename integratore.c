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
  double x0, y0, z0;		// 
  FILE *output;			    // def un puntatore per il file out
  FILE *input;
  FILE *periodo;

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

/****Creo un output.dat e scrivo l'ordine dei dati****/
  output = fopen("output.dat", "w"); //iniz il file
  fprintf(output, "#t x y z\n");

/****Apro un altro file per i periodi****/
  periodo = fopen("perio.dat", "w");
  fprintf(periodo, "#periodo ");

  //x0, y0 e z0
  x0 = p.x;
  y0 = p.y;
  z0 = p.z;

  fprintf(output, "%.8lf %.8lf %.8lf %.8lf\n", 0., x0, y0, z0);

	steps = (tmax/c.dt); //assegno a steps il numero di passi di integrazione


/****Zona principale del codice****/
	for (i=0; i<steps; i++) {  
    
    pold = strcopy(p);  //ricopio il vettore dentro uno nuovo per poter controllare il periodo

    p = RK4(p, c, (double)(i + 1) * c.dt);
    fprintf(output, "%.4lf %.8lf %.8lf %.8lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);


  /****Zona del controllo periodo****/
    //if (p.x <= pold.x && p.y - y0 <= y0*0.001 && p.z - z0 <= z0*0.001) { //se periodica dovrebbe tornare al punto di partenza
    if (p.x >= x0 && pold.x <= x0 && p.y >= y0 && pold.y <= y0 && p.z >= z0 && pold.z <= z0) {
      dperiodo[k] = c.dt*((double)(i+1)); //temporaneamente metto i tempi in un array
      k++;  //conto quanti ne salvo
      // /****REALLOC ZONE - in fase di costruzione****/
      // if (k+1 == NUM_PERIODI * (j+1)) { //k+1, prima di finire completamente la memoria ne assegno della nuova
      //   dperiodo = realloc(dperiodo, 2 * NUM_PERIODI * sizeof(double)); //raddoppio lo spazio riservato a dperiodo
      //   j++; //j mi aiuta a tenere ordinata la memoria
      // }
    }

	}


/****CALCOLO E SALVATAGGIO DEL PERIODO****/
  fprintf(periodo, "%d\n", k);
  
  for (i=0; i<k-1; i++) {
    fprintf(periodo, "%lf\n", dperiodo[i+1]-dperiodo[i]);
  }


/****CONCLUSIONI****/

  fclose(input);
  free(dperiodo); //heap libero
  fclose(output);
  fclose(periodo);

	exit(0);
}



struct vector RK4(vector n, pars c, double t) { //chiamo n il vettore così n mi richiama alla mente a che passo di int sono
  //i passi xi, yi e zi sono le velocità
  
  double x1, x2, x3, x4;  //definisco i passi di RK4 temporanee per x
  double dt = c.dt;

  x1 = f(n.x, t, n, c) * dt;
  x2 = f(n.x + 0.5 * x1, t + 0.5 * dt, n, c) * dt;    //repetita iuvant: i passi vengono fatti per tutte e tre le dimensioni
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

  n.z += (z1 + 2.*z2 + 2.*z3 + z4)/6.;

  return n;
}



/****LE FUNZIONI f, g e h del problema(in ordine x, y, z)****/

double f(double x, double t, vector n, pars c) {
  return x * ( c.a * (1 - n.x) + 0.5 * (1 - n.y) + c.b * (1 - n.z) );
}

double g(double y, double t, vector n, pars c) {
  return y * ( -0.5 * (1 - n.x) + c.b * (n.y - n.z) );
}

double h(double z, double t, vector n, pars c) {
  return z * ( c.rho * (1 - n.x) + 0.1 * (2. - n.y - n.z) );
}



/****Utili****/
struct vector strcopy(vector a) {
  
  vector new;
  
  new.x = a.x;
  new.y = a.y;
  new.z = a.z;

  return new;
}