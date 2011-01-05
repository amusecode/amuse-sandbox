#include <stdio.h>
#include "kepler.h"

int main (int argc, char *argv[]) {

  double t;
  double x=0.0;
  double y=0.0;
  double z=0.0;
  char X[600];
  char Y[600];
  char Z[600];
  char m[300];

  initialize(100);
  set_position(1.0, 0.1, -0.1);
  set_velocity(-0.1, 0.1, -0.2);
  set_mu(1.0);
  for (t = 0.01; t < 20; t+=0.005) {
    if (evolve_d(t)==-1) {
      printf("WARNING: Newton root finding failed @t=%2.3e\n", t);
    }
    else {
      get_position_s(X, Y, Z, 10);
      sprintf(m, "%s %s %s\n", X, Y,Z);
      printf(m);
      get_position(&x, &y, &z);
    }
  }
  return 0;
}
