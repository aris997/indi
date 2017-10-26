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

int asintotico(vector, vector, pars);

int main(int argv, char *argc[]) {

	int tmax, k=0, j=0;
	long int i, steps;
	double x0, y0, z0, C, Co, e, per;
	double dperiodo[256];

	FILE *output;
	FILE *input;
	FILE *periodo;
	FILE *cost;
	FILE *piripicchio;

	vector p, pold;
	pars c;

	input = fopen("input.dat", "r");
	fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &p.x, &p.y, &p.z, &c.dt, &tmax);

	output = fopen("output.dat", "w");
	fprintf(output, "#t x y z\n");

	periodo = fopen("perio.dat", "w");
	fprintf(periodo, "#periodo\n");

	cost = fopen("costante.dat", "w");
	fprintf(cost, "#Costante\n");

	piripicchio = fopen("piripicchio.dat", "w");

	x0 = p.x;
	y0 = p.y;
	z0 = p.z;

	e = c.dt/2.;

	fprintf(output, "%.8lf %.8lf %.8lf %.8lf\n", 0., x0, y0, z0);

	Co = log(x0) - x0 + log(y0) - y0;

	steps = (tmax/c.dt);

	for (i=0; i<steps; i++) {

		pold = strcopy(p);

		p = RK4(p, c);
		fprintf(output, "%.8lf %.16lf %.16lf %.16lf\n", c.dt*((double)(i+1)), p.x, p.y, p.z);

		C = log(p.x) - p.x + log(p.y) - p.y + p.z;
		fprintf(cost, "%lf %.16lf\n", (double)i*c.dt, C);

		if ( (double)i*c.dt == 10.) {
			printf("%.52lf ", C/Co - 1);
			fprintf(piripicchio, "%.48lf \n", C/Co - 1);
		}

		// if () {
		// 	dperiodo[k] = (double)i*c.dt;
		// 	printf("%lf %lf %lf\n", pold.x, p.x, (double)i*c.dt);
		// 	k++;
		// }
	}
	

	for (i=0; i<k-1; i++) {
		per = dperiodo[i+1]-dperiodo[i];
		if (per > c.dt*100.) fprintf(periodo, "%lf\n", dperiodo[i+1]-dperiodo[i]); //eseguo una pulizia dati
	}





	FILE *rho_variable;
	char filename[30];

	steps = (long int)(tmax/c.dt);


	for (j=0; j<2; j++) {

			fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d\n", &c.a, &c.b, &c.rho, &x0, &y0, &z0, &c.dt, &tmax);
			sprintf(filename, "rho%d.dat", (int)(c.rho*100.));

			//printf("%d ", (int)(c.rho*100.));

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

	FILE *drastico;
	FILE *gpscript;

	gpscript = fopen("mscript.gp", "w");
	fprintf(gpscript, "set term x11 persist\nset multiplot\nunset autoscale\nset yrange [0.5:8]\nset xrange [200:400]\nunset key\nplot");

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
	int width = 3;
	vector S;

	steps = (long int)(tmax/c.dt);
	//printf("%ld\n", steps);


	for (j=0; j<=100; j++)  { //100 passi per arrivare da 1.2 a 1.5 rho


		c.rho = 1.2;
		c.rho += 0.003*(double)j; //incremento dell'1% dell'intervallo

		sprintf(filename, "cambiamento_drast_rho%0*d.dat", width, j);
		fprintf(gpscript, " '%s' every ::5000::8000 u 1:4 w l\nreplot", filename);
		drastico = fopen(filename, "w");

		sprintf(filename, "periodo_drast_rho%0*d.dat", width, j);
		//completare la scrittura dello script gp
		periodo = fopen(filename, "w");
		//printf("%lf \n", c.rho);

		p.x = x0;
		p.y = y0;
		p.z = z0;

		for (i = 0; i<=steps; i++) { //steps dovrebbe essere 8000

			pold = strcopy(p);

			p = RK4(p, c);

			//if (i%200 == 0) printf("%.8lf %.16lf %.16lf %.16lf ", (double)i*c.dt, p.x, p.y, p.z);
			if (i == 4000) S = strcopy(p);
			if (i > 4000) {
				fprintf(drastico, "%.8lf %.16lf %.16lf %.16lf\n", (double)i*c.dt, p.x, p.y, p.z);
				if (e >= fabs(p.x - S.x) && e >= fabs(p.y - S.y) && e >= fabs(p.z - S.z)) {
					fprintf(periodo, "%s\n", );
				}
			}

			
			//As += asintotico(pold, p, c);
			//if (As == 10000) { i = steps; printf("%lf esiste asintoto\n", c.rho); }
		}

		fclose(drastico);

	}


	fclose(periodo);
	fclose(cost);
	fclose(input);
	fclose(output);
	fclose(piripicchio);

	exit(0);
}


int asintotico(vector pold, vector p, pars c) {

	return 1;
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