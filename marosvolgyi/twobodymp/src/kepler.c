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

void stumpff(mpfr_t s, mpfr_t  c0, mpfr_t c1, mpfr_t c2, mpfr_t c3) {
  //min = -x, inv = 1/x, neg = 1-x

  mpfr_t sqrt_s;
  mpfr_t my_cos, my_sin, my_divsin;
  mpfr_t min_s, inv_s, neg_s;
  int condition;

  mpfr_inits2(PR, sqrt_s, my_cos, my_sin, my_divsin, 
	      min_s, inv_s, neg_s,
	      (mpfr_ptr)0);
  condition = mpfr_cmp_si(s, 0);

  if (condition>0) {
    mpfr_sqrt(sqrt_s, s, GMP_RNDN);
    mpfr_cos(my_cos, sqrt_s, GMP_RNDN);
    mpfr_sin(my_sin, sqrt_s, GMP_RNDN);
    mpfr_div(my_divsin, my_sin, sqrt_s, GMP_RNDN);
    mpfr_set(c0, my_cos, GMP_RNDN);
    mpfr_set(c1, my_divsin, GMP_RNDN);
  }
  else if (condition<0){
    mpfr_neg(min_s, s, GMP_RNDN);
    mpfr_sqrt(sqrt_s, min_s, GMP_RNDN);
    mpfr_cosh(my_cos, sqrt_s, GMP_RNDN);
    mpfr_sinh(my_sin, sqrt_s, GMP_RNDN);
    mpfr_div(my_divsin, my_sin, sqrt_s, GMP_RNDN);
    mpfr_set(c0, my_cos, GMP_RNDN);
    mpfr_set(c1, my_divsin, GMP_RNDN);
  }
  else printf("Error in stumpff s = 0\n");

  mpfr_d_div(inv_s, 1.0, s, GMP_RNDN);
  mpfr_si_sub(neg_s, 1, c0, GMP_RNDN);
  mpfr_mul(c2, inv_s, neg_s, GMP_RNDN);
  mpfr_si_sub(neg_s, 1, c1, GMP_RNDN);
  mpfr_mul(c3, inv_s, neg_s, GMP_RNDN);

  mpfr_clears(sqrt_s, my_cos, my_sin, my_divsin, min_s, inv_s, neg_s, (mpfr_ptr)0);//checked!
}

int initial_guess_for_s (mpfr_t dt, 
			 mpfr_t r0, mpfr_t u, mpfr_t alpha,
			 mpfr_t s) {
  mpfr_t A, En, Ec, Es, E, X, Y, Dm;
  mpfr_t dtr0_abs, dt_over_r0;
  mpfr_t foo, bar, dttemp;
  mpfr_t pi, dpi;//dpi is double pi = 2*pi
  int crit;
  short sigma;

  mpfr_inits2(PR, A, En, Ec, Es, E, X, Y, Dm, 
	      dtr0_abs,
	      dt_over_r0,
	      foo,
	      bar,
	      dttemp,
	      pi,
	      dpi,
	      (mpfr_ptr)0);

  mpfr_div(dt_over_r0, dt, r0, GMP_RNDN);
  mpfr_abs(dtr0_abs, dt_over_r0, GMP_RNDN);
  crit = mpfr_cmp_d(dtr0_abs, 0.2);
  if (crit < 0) {
    mpfr_mul(foo, dt, dt, GMP_RNDN);
    mpfr_mul(foo, foo, u, GMP_RNDN);
    mpfr_div_d(foo, foo, 2.0, GMP_RNDN);
    mpfr_div(foo, foo, r0, GMP_RNDN);
    mpfr_div(foo, foo, r0, GMP_RNDN);
    mpfr_div(foo, foo, r0, GMP_RNDN);
    mpfr_sub(s, dt_over_r0, foo, GMP_RNDN); 
  }
  else if (alpha>0) {
    //elliptic motion initial guess
    mpfr_div(A, mu, alpha, GMP_RNDN);

    mpfr_div(foo, mu, A, GMP_RNDN);
    mpfr_div(foo, foo, A, GMP_RNDN);
    mpfr_div(foo, foo, A, GMP_RNDN);
    mpfr_sqrt(En, foo, GMP_RNDN);

    mpfr_div(foo, r0, A, GMP_RNDN);
    mpfr_d_sub(Ec, 1.0, foo, GMP_RNDN);

    mpfr_div(foo, u, En, GMP_RNDN);
    mpfr_div(foo, foo, A, GMP_RNDN);
    mpfr_div(Es, foo, A, GMP_RNDN);

    mpfr_pow_si(foo, Ec, 2, GMP_RNDN);
    mpfr_pow_si(bar, Es, 2, GMP_RNDN);
    mpfr_add(foo, foo, bar, GMP_RNDN);
    mpfr_sqrt(E, foo, GMP_RNDN);
    
    mpfr_const_pi(pi, GMP_RNDN);
    mpfr_mul_d(dpi, pi, 2.0, GMP_RNDN);
    mpfr_div(foo, dpi, En, GMP_RNDN);
 
    mpfr_mul(bar, En, dt, GMP_RNDN);
    mpfr_div(bar, bar, dpi, GMP_RNDN);
    long flooor = mpfr_get_si(bar, GMP_RNDN);
    mpfr_mul_si(foo, foo, flooor, GMP_RNDN);
    mpfr_sub(dttemp, dt, foo, GMP_RNDN);

    mpfr_mul(foo, En, dttemp, GMP_RNDN);
    mpfr_sub(Y, foo, Es, GMP_RNDN);

    mpfr_sin(foo, Y, GMP_RNDN);
    mpfr_mul(foo, Ec, foo, GMP_RNDN);
    mpfr_cos(bar, Y, GMP_RNDN);
    mpfr_mul(bar, Es, bar, GMP_RNDN);
    mpfr_add(foo, foo, bar, GMP_RNDN);
    sigma = mpfr_sgn(foo);
    mpfr_mul_d(foo, E, 1.0*sigma*0.85, GMP_RNDN);
    mpfr_add(X, foo, Y, GMP_RNDN);

    mpfr_sqrt(foo, alpha, GMP_RNDN);
    mpfr_div(s, X, foo, GMP_RNDN);
  }
  else {
    //hyperbolic motion
    mpfr_div(A, mu, alpha, GMP_RNDN);
    mpfr_div(foo, mu, A, GMP_RNDN);
    mpfr_div(foo, foo, A, GMP_RNDN);
    mpfr_div(foo, foo, A, GMP_RNDN);
    mpfr_neg(foo, foo, GMP_RNDN);
    mpfr_sqrt(En, foo, GMP_RNDN);
    
    mpfr_div(foo, r0, A, GMP_RNDN);
    mpfr_si_sub(Ec, 1, foo, GMP_RNDN);

    mpfr_mul(foo, A, mu, GMP_RNDN);
    mpfr_neg(foo, foo, GMP_RNDN);
    mpfr_sqrt(foo, foo, GMP_RNDN);
    mpfr_div(Es, u, foo, GMP_RNDN); 

    mpfr_mul(foo, Ec, Ec, GMP_RNDN);
    mpfr_mul(bar, Es, Es, GMP_RNDN);
    mpfr_sub(foo, foo, bar, GMP_RNDN);
    mpfr_sqrt(E, foo, GMP_RNDN);

    mpfr_mul(Dm, En, dt, GMP_RNDN);

    if (mpfr_cmp_si(Dm, 0)<0) {
      mpfr_mul_d(foo, Dm, -2.0, GMP_RNDN);
      mpfr_mul_d(bar, E, 1.8, GMP_RNDN);
      mpfr_add(foo, foo, bar, GMP_RNDN);
      mpfr_sub(bar, Ec, Es, GMP_RNDN);
      mpfr_div(foo, foo, bar, GMP_RNDN);
      mpfr_log(foo, foo, GMP_RNDN);
      mpfr_neg(foo, foo, GMP_RNDN);
      mpfr_neg(bar, alpha, GMP_RNDN);
      mpfr_sqrt(bar, bar, GMP_RNDN);
      mpfr_div(s, foo, bar, GMP_RNDN);
    }
    else {
      mpfr_mul_d(foo, Dm, 2.0, GMP_RNDN);
      mpfr_mul_d(bar, E, 1.8, GMP_RNDN);
      mpfr_add(foo, foo, bar, GMP_RNDN);
      mpfr_add(bar, Ec, Es, GMP_RNDN);
      mpfr_div(foo, foo, bar, GMP_RNDN);
      mpfr_log(foo, foo, GMP_RNDN);
      mpfr_neg(bar, alpha, GMP_RNDN);
      mpfr_sqrt(bar, bar, GMP_RNDN);
      mpfr_div(s, foo, bar, GMP_RNDN);
    }
  }

  mpfr_clears(A, En, Ec, Es, E, X, Y, Dm, 
	      dtr0_abs, dt_over_r0, foo, bar, dttemp, pi, dpi,
	      (mpfr_ptr)0);//checked!

  return 0;
}

int kepler_solve (mpfr_t dt, 
		  mpfr_t F, mpfr_t G, mpfr_t Fdot, mpfr_t Gdot) {
  mpfr_t f, fp;
  mpfr_t ds, s;
  mpfr_t r0, v0s, u, alpha;
  mpfr_t c0, c1, c2, c3;
  mpfr_t ssalpha, foo, bar, dummy;
  int i=0;

  mpfr_inits2(PR, f, fp, ds, s,
	      r0, v0s, u, alpha,
	      c0, c1, c2, c3,
	      ssalpha, foo, bar, dummy,
	      (mpfr_ptr)0);

  mpfr_mul(foo, x,x, GMP_RNDN);
  mpfr_mul(bar, y,y, GMP_RNDN);
  mpfr_add(foo, foo, bar, GMP_RNDN);
  mpfr_mul(bar, z, z, GMP_RNDN);
  mpfr_add(foo, foo, bar, GMP_RNDN);
  mpfr_sqrt(r0, foo, GMP_RNDN);

  mpfr_mul(foo, vx, vx, GMP_RNDN);
  mpfr_mul(bar, vy, vy, GMP_RNDN);
  mpfr_add(foo, foo, bar, GMP_RNDN);
  mpfr_mul(bar, vz, vz, GMP_RNDN);
  mpfr_add(v0s, foo, bar, GMP_RNDN);

  mpfr_mul(foo, x, vx, GMP_RNDN);
  mpfr_mul(bar, y, vy, GMP_RNDN);
  mpfr_add(foo, foo, bar, GMP_RNDN);
  mpfr_mul(bar, z, vz, GMP_RNDN);
  mpfr_add(u, foo, bar, GMP_RNDN);

  mpfr_mul_d(foo, mu, 2.0, GMP_RNDN);
  mpfr_div(foo, foo, r0, GMP_RNDN);
  mpfr_sub(alpha, foo, v0s, GMP_RNDN);

  initial_guess_for_s (dt, r0, u, alpha, s);
  mpfr_set_d(ds, 1.0, GMP_RNDN);
  mpfr_abs(foo, ds, GMP_RNDN);
  
  while (mpfr_cmp_d(foo,1.0e-24 ) > 0) {
    mpfr_mul(foo, s, s, GMP_RNDN);
    mpfr_mul(ssalpha, foo, alpha, GMP_RNDN);
    stumpff(ssalpha, c0, c1, c2, c3);
    mpfr_mul(c1, s, c1, GMP_RNDN);
    mpfr_mul(c2, s, c2, GMP_RNDN);
    mpfr_mul(c2, s, c2, GMP_RNDN);
    mpfr_mul(c3, s, c3, GMP_RNDN);
    mpfr_mul(c3, s, c3, GMP_RNDN);
    mpfr_mul(c3, s, c3, GMP_RNDN);

    mpfr_mul(foo, r0, c1, GMP_RNDN);
    mpfr_mul(bar, u, c2, GMP_RNDN);
    mpfr_add(foo, foo, bar, GMP_RNDN);
    mpfr_mul(bar, mu, c3, GMP_RNDN);
    mpfr_add(foo, foo, bar, GMP_RNDN);
    mpfr_sub(f, foo, dt, GMP_RNDN);

    mpfr_mul(foo, r0, c0, GMP_RNDN);
    mpfr_mul(bar, u, c1, GMP_RNDN);
    mpfr_add(foo, foo, bar, GMP_RNDN);
    mpfr_mul(bar, mu, c2, GMP_RNDN);
    mpfr_add(fp, foo, bar, GMP_RNDN);

    mpfr_div(ds, f, fp, GMP_RNDN);
    mpfr_sub(s, s, ds, GMP_RNDN);

    if (i++>50) {
      printf("Convergence error in Newton method\n");
      return -1; 
    }
    mpfr_abs(foo, ds, GMP_RNDN);
  
  }

  mpfr_div(foo, mu, r0, GMP_RNDN);
  mpfr_mul(foo, foo, c2, GMP_RNDN);
  mpfr_d_sub(F, 1.0, foo, GMP_RNDN);

  mpfr_mul(foo, mu, c3, GMP_RNDN);
  mpfr_sub(G, dt, foo, GMP_RNDN);
  
  mpfr_div(foo, mu, fp, GMP_RNDN);
  mpfr_div(foo, foo, r0, GMP_RNDN);
  mpfr_mul(foo, foo, c1, GMP_RNDN);

  mpfr_neg(Fdot, foo, GMP_RNDN);

  mpfr_div(foo, mu, fp, GMP_RNDN);
  mpfr_mul(foo, foo, c2, GMP_RNDN);
  mpfr_si_sub(Gdot, 1.0, foo, GMP_RNDN);

  mpfr_clears(f, fp, ds, s,
	      r0, v0s, u, alpha,
	      c0, c1, c2, c3,
	      ssalpha, foo, bar, dummy,
	      (mpfr_ptr)0);

  return 0;
}

int evolve_d(double time_new_) {
  mpfr_t time_new;
  int result;

  mpfr_init2(time_new, PR);
  mpfr_set_d(time_new, time_new_, GMP_RNDN);

  result = evolve(time_new);

  mpfr_clear(time_new);
  return result;
}
  
int evolve (mpfr_t time_new) {
  mpfr_t x_new, y_new, z_new;
  mpfr_t vx_new, vy_new, vz_new;
  mpfr_t F, G, Fdot, Gdot;
  mpfr_t dt;
  mpfr_t foo, bar;
  
  mpfr_inits2(PR, 
	      x_new, y_new, z_new,
	      vx_new, vy_new, vz_new,
	      F, G, Fdot, Gdot,
	      dt,
	      foo,
	      bar,
	      (mpfr_ptr)0);

  mpfr_sub(dt, time_new, time, GMP_RNDN);

  if (kepler_solve(dt, F, G, Fdot, Gdot) == 0) { 
    mpfr_mul(foo, x, F, GMP_RNDN);
    mpfr_mul(bar, vx, G, GMP_RNDN);
    mpfr_add(x_new, foo, bar, GMP_RNDN);

    mpfr_mul(foo, y, F, GMP_RNDN);
    mpfr_mul(bar, vy, G, GMP_RNDN);
    mpfr_add(y_new, foo, bar, GMP_RNDN);

    mpfr_mul(foo, z, F, GMP_RNDN);
    mpfr_mul(bar, vz, G, GMP_RNDN);
    mpfr_add(z_new, foo, bar, GMP_RNDN);

    mpfr_mul(foo, x, Fdot, GMP_RNDN);
    mpfr_mul(bar, vx, Gdot, GMP_RNDN);
    mpfr_add(vx_new, foo, bar, GMP_RNDN);

    mpfr_mul(foo, y, Fdot, GMP_RNDN);
    mpfr_mul(bar, vy, Gdot, GMP_RNDN);
    mpfr_add(vy_new, foo, bar, GMP_RNDN);

    mpfr_mul(foo, z, Fdot, GMP_RNDN);
    mpfr_mul(bar, vz, Gdot, GMP_RNDN);
    mpfr_add(vz_new, foo, bar, GMP_RNDN);

    mpfr_set(x, x_new, GMP_RNDN);
    mpfr_set(y, y_new, GMP_RNDN);
    mpfr_set(z, z_new, GMP_RNDN);
    mpfr_set(vx, vx_new, GMP_RNDN);
    mpfr_set(vy, vy_new, GMP_RNDN);
    mpfr_set(vz, vz_new, GMP_RNDN);

    mpfr_set(time, time_new, GMP_RNDN);
    return 0;
  }
  else return -1;

  mpfr_clears(x_new, y_new, z_new,
	      vx_new, vy_new, vz_new,
	      F, G, Fdot, Gdot,
	      dt,
	      foo, bar,
	      (mpfr_ptr)0);
  mpfr_free_cache ();
}

int initialize(int precision) {
  if (precision<MPFR_PREC_MIN) {
    precision = MPFR_PREC_MIN;
    fprintf(stderr, "Warning: set precision to minimal precision\n");
  }
  if (precision>MPFR_PREC_MAX) {
    precision = MPFR_PREC_MAX;
    fprintf(stderr, "Warning: set precision to maximal precision\n");
  }

  mpfr_set_default_prec (precision);
  PR = precision;

  mpfr_inits2(PR, time,
	      x, y, z, vx, vy, vz,
	      mu, (mpfr_ptr)0);

  mpfr_set_d(time, 0.0, GMP_RNDN);
  mpfr_set_d(x, 0.0, GMP_RNDN);
  mpfr_set_d(y, 0.0, GMP_RNDN);
  mpfr_set_d(z, 0.0, GMP_RNDN);
  mpfr_set_d(vx, 0.0, GMP_RNDN);
  mpfr_set_d(vy, 0.0, GMP_RNDN);
  mpfr_set_d(vz, 0.0, GMP_RNDN);
  mpfr_set_d(mu, 0.0, GMP_RNDN);
  if (strcmp (mpfr_get_version (), MPFR_VERSION_STRING))
    fprintf (stderr, "Warning: header and library do not match\n");
  
  id = 0;

  return 0;
}

int set_mu(double mu_) {
  mpfr_set_d(mu, mu_, GMP_RNDN);
  return 0;
}

int set_pos(double x_, double y_, double z_) {
  mpfr_set_d(x, x_, GMP_RNDN);
  mpfr_set_d(y, y_, GMP_RNDN);
  mpfr_set_d(z, z_, GMP_RNDN);
  return 0;
}

int get_pos(double *x_, double *y_, double *z_) {
  *x_ = mpfr_get_d(x, GMP_RNDN);
  *y_ = mpfr_get_d(y, GMP_RNDN);
  *z_ = mpfr_get_d(z, GMP_RNDN);
  return 0;
}

int get_position_s(char *X, char *Y, char *Z, int pr) {
  char format[10];
  sprintf(format, "%%.%dRf", pr);
  mpfr_sprintf (X, format, x);
  mpfr_sprintf (Y, format, y);
  mpfr_sprintf (Z, format, z);
  return 0;
}

int set_vel(double vx_, double vy_, double vz_) {
  mpfr_set_d(vx, vx_, GMP_RNDN);
  mpfr_set_d(vy, vy_, GMP_RNDN);
  mpfr_set_d(vz, vz_, GMP_RNDN);
  return 0;
}

int get_vel(double *vx_, double *vy_, double *vz_) {
  *vx_ = mpfr_get_d(vx, GMP_RNDN);
  *vy_ = mpfr_get_d(vy, GMP_RNDN);
  *vz_ = mpfr_get_d(vz, GMP_RNDN);
  return 0;
}

void test_precision() {
  mpfr_t pi;
  mpfr_init2(pi, 512);
  mpfr_const_pi(pi, GMP_RNDN);
  mpfr_printf("pi = %.512RNf\n", pi);
  fflush(stdout);
  mpfr_clear(pi);
}
