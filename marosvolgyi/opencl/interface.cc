/* 
   AMUSE interface to opencl nbody code
   ====================================
   
 */

#include "worker_code.h"

#ifdef __PLOT
#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
}
#endif
#endif

int mainentry(int argc, char *argv[]);

int initialization() 
{
  char *params[] = {"nbody", "-x", " 4000"};
  mainentry(3, params);
  return 0;
}

int new_particle(int *id_, double mass_, double radius, 
		 double x_, double y_, double z_,
                 double vx_, double vy_, double vz_) 
{
  return 0;
}

int get_number_of_particles(int *number_of_particles) 
{
  return 0;
}

int commit_particles() 
{
  return -1;
}

int evolve_system(double new_time) 
{
  return 0;
}

