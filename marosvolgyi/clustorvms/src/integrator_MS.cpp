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

#include "/home/boekholt/Packages/mpfr_c++/mpreal.h"
using namespace mpfr;

#include "STAR.h"
#include "FORCE.h"
#include "DYNAMICS.h"
#include "CLOCK.h"
#include "CLUSTER.h"
#include "BS_INTEGRATOR.h"

int main()
{
  ////////////////////////////////////////////////////////
  // Read in configuration from input file
  ////////////////////////////////////////////////////////
  string file_in, file_out;
  int n_max, k_max, numBits;
  mpreal t_sim, dt_min, dt_max, dv_max, dt_print, epsilon, softening;  

  string dummy;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> file_in;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> file_out;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> t_sim;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> dt_min;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> dt_max;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> dv_max;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> dt_print;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> numBits;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> epsilon;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> n_max;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> k_max;
  for(int i=0; i<3; i++) cin >> dummy;
  cin >> softening;    

  ////////////////////////////////////////////////////////
  // Set mpreal memory
  ////////////////////////////////////////////////////////
  mpreal::set_default_prec(numBits);

  ////////////////////////////////////////////////////////
  // Configure Objects
  ////////////////////////////////////////////////////////
  FORCE force(softening);
  CLUSTER cluster(file_in, force);
  BS_INTEGRATOR bs(epsilon, n_max, k_max);
  CLOCK clock(cluster.get_t(), t_sim, dt_min, dt_max, dv_max, dt_print);

  ////////////////////////////////////////////////////////
  // Run Simulation
  ////////////////////////////////////////////////////////
  clock.start_timer();

  ofstream odata;
  odata.precision(numBits/4);
  odata.open( file_out.c_str() );
  if( !odata )
  {
    cerr << "Could not open " << file_out << "!" << endl;
    exit(1);
  }
  else
  {
    cluster.print(odata);
    cerr << endl;
    cerr << clock.get_progress() << "%" << endl;

    while( !clock.alarm() )
    {
      cluster.calc_a();
      clock.calc_dt( cluster.get_a_max() );

      bs.integrate(cluster, clock.get_dt());

      if( !bs.converged() )
      {
        cerr << "No Convergence Reached, Simulation Aborted!" << endl;
        clock.abort();
      }
      else
      {
        clock.tick();
        cluster.set_t( clock.get_t() );
        if( clock.to_print() )
        {
	  CLUSTER cl_exp = cluster;
          cl_exp.calc_a();
	  mpreal dt_exp = clock.get_t_print() - clock.get_t();
	  bs.integrate(cl_exp, dt_exp);
	  cl_exp.set_t( clock.get_t_print() );
          cl_exp.print(odata);
        }
      }
      cerr << clock.get_progress() << "%" << endl;
    }
  }
  odata.close();  

  clock.stop_timer();

  ////////////////////////////////////////////////////////
  // Finish Program
  ////////////////////////////////////////////////////////
  cerr << endl;
  cerr << "Runtime: " << clock.get_timer() << endl;
  cerr << endl;  

  return 0;
}
