#include <stdio.h>

int evolve(double t);
int set_position(double r[3]);
int get_position(double r[3]);
int set_velocity(double v[3]);

int main (int argc, char *argv[]) {
  double t;
  double r[3];
  double v[3];
  FILE *datafile;
  
  datafile = fopen("kepler.dat","w");
  printf("Writing result to file 'kepler.dat'\n");
  r[0] = 1.0; r[1] = 0.1; r[2] = -0.1;
  v[0] = -0.1; v[1] = 2.0; v[2] = -0.2;
  r[0] = 1.0; r[1] = 0.0; r[2] = 0.0;
  v[0] = 0.0; v[1] = .01; v[2] = 0.0;

  set_position(r);
  set_velocity(v);
  set_mu(1.0);
  for (t = 0.01; t<5.0; t+=0.05) {
    if (evolve(t)==-1) {
      printf("WARNING: Newton root finding failed @t=%2.3e\n", t);
    }
    else{
      get_position(r);
      fprintf(datafile, "%2.3e, %2.3e, %2.3e\n", r[0], r[1], r[2]);
    }
  }
  fflush(datafile);
  fclose(datafile);
  return 0;
}
