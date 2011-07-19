//interface to:

//////////////////////////////////////////////////////////////
// integrator_MS = N-Body code using Mpreals and Single core
//
// Input  = fixed format file with parameters
// Output = textfile with newly generated states, with path defined in Input
//
// Tjarda Boekholt
// July 2011
//////////////////////////////////////////////////////////////

#include <iostream>
using namespace std;

#include <fstream>
#include <string>

#include <sys/time.h>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "./Packages/mpfr_c++/mpreal.h"
using namespace mpfr;

/*
#include "STAR.h"
#include "FORCE.h"
#include "DYNAMICS.h"
#include "CLOCK.h"
#include "CLUSTER.h"
#include "BS_INTEGRATOR.h"
*/
string file_in, file_out;
int n_max, k_max, numBits;
mpreal t_sim, dt_min, dt_max, dv_max, dt_print, epsilon, softening;  

int echo(int input){
    return input;
}


