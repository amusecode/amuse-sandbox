/* 
Multiprecision two-body problem solver

code rewrite 2 c and mp by Marcell Marosvolgyi 
marosvolgyi@strw.leidenuniv.nl

code based on example code in:
Fundamentals of Celestial Mechanics, J.M.A. Danby 2nd Edition
*/

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

int PR;

mpfr_t x, y, z;
mpfr_t vx, vy, vz;
mpfr_t time ;
mpfr_t mu;
int id;

#ifdef __cplusplus
extern "C" {
#endif

void stumpff(mpfr_t s, mpfr_t  c0, mpfr_t c1, mpfr_t c2, mpfr_t c3); 

int initial_guess_for_s (mpfr_t dt, 
			 mpfr_t r0, mpfr_t u, mpfr_t alpha,
			 mpfr_t s);

int kepler_solve (mpfr_t dt, 
		  mpfr_t F, mpfr_t G, mpfr_t Fdot, mpfr_t Gdot);

int evolve_d(double time_new_);

int evolve (mpfr_t time_new);

int initialize(int precision);

int set_mu(double mu_);

int set_pos(double x_, double y_, double z_);

int get_pos(double *x_, double *y_, double *z_); 

int get_position_s(char *X, char *Y, char *Z, int pr);

int get_velocity_s(char *VX, char *VY, char *VZ, int pr);

int set_vel(double vx_, double vy_, double vz_);

int get_vel(double *vx_, double *vy_, double *vz_);

void test_precision();

#ifdef __cplusplus
}
#endif
