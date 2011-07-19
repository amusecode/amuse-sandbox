extern int echo(int input);
#include <iostream>
using namespace std;

#include <fstream>
#include <string>

#include <sys/time.h>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "./src/Packages/mpfr_c++/mpreal.h"
using namespace mpfr;


#include "./src/STAR.h"
#include "./src/FORCE.h"
#include "./src/DYNAMICS.h"
#include "./src/CLOCK.h"
#include "./src/CLUSTER.h"
#include "./src/BS_INTEGRATOR.h"

/*
 * Interface code
 */
 
int number_of_stars = 0;
CLUSTER cluster();

int echo_int(int input, int * output){
  *output = echo(input);
    return 0;
}

int set_numBits(int number_of_bits)
{
  mpreal::set_default_prec(number_of_bits);
  return 0;
}

int set_softening(double softening)
{
  FORCE force(softening);
  return 0;
}

int new_particle(int *id, 
		 double mass, double radius,
		 double x, double y, double z, 
		 double vx, double vy, double vz)
{
  //STAR st(mass, x, y, z, vx, vy, vz);
  //....star.push_back(st);
  //*id = ++number_of_stars;
  return 0;
}
