/* 
Kepler integrator, two-body problem

code based on:
Fundamentals of Celestial Mechanics, J.M.A. Danby 2nd Edition
*/

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

#define PR 10 //precision

mpfr_t x, y, z;
mpfr_t vx, vy, vz;
mpfr_t time = 0.0;
mpfr_t mu;

int sign(mpfr_t x) {
  mpfr_t xabs;
  int result;

  mpfr_init2(xabs, PR);
  mpfr_abs(xabs,x,GMP_RNDN);
  if (mpfr_cmp_si(xabs,0)<0) result =  -1;
  else result = 1;
  mpfr_clear(xabs);
  return result;
}

void stumpff(mpfr_t s, mpfr_t  c0, mpfr_t c1, mpfr_t c2, mpfr_t c3) {
  mpfr_t sqrt_s;
  mpfr_t my_cos, my_sin, my_divsin;
  mpfr_t min_s, inv_s;
  int condition;

  mpfr_inits2(PR, sqrt_s, my_cos, my_sin, my_divsin, min_s, inv_s, (mpfr_ptr)0);
  condition = mpfr_comp_si(s, 0);

  if (condition>0) {
    mpfr_sqrt(sqrt_s, s, GMP_RNDN);
    mpfr_cos(my_cos, sqrt_s, GMP_RNDN);
    mpfr_sin(my_sin, sqrt_s, GMP_RNDN);
    mpfr_div(my_divsin, my_sin, sqrt_s, GMP_RNDN);
    c0 = my_cos;
    c1 = my_divsin;
  }
  else if (condition<0){
    mpfr_neg(min_s, s, GMP_RNDN);
    mpfr_sqrt(sqrt_s, min_s, GMP_RNDN);
    mpfr_cosh(my_cos, sqrt_s, GMP_RNDN);
    mpfr_sinh(my_sin, sqrt_s, GMP_RNDN);
    mpfr_div(my_divsin, my_sin, sqrt_s, GMP_RNDN);
    c0 = my_cos;
    c1 = my_divsin;
  }
  else printf("Error in stumpff s = 0\n");

  mpfr_d_div(inv_s, 1.0, s, GMP_RNDN);
  //cont here...

  *c2 = inv_s * (1.0 - *c0);
  *c3 = inv_s * (1.0 - *c1);

  mpfr_clears(sqrt_s, my_cos, my_sin, my_divsin, min_s, inv_s, (mpfr_ptr)0);
}

double initial_guess_for_s (mpfr_t dt, 
			    mpfr_t r0, mpfr_t u, mpfr_t alpha) {
  double A, En, Ec, Es, E, X, Y, Dm, sigma;
  double s;

  if (mpfr_abs(dt/r0, dt/r0, GMP_RNDN) <=.2) {
    s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0);
  }
  else if (alpha>0) {
    //elliptic motion initial guess
    A = mu/alpha;
    En = pow(mu/A/A/A, 0.5);
    Ec = 1.0 - r0/A;
    Es = u/En/A/A;
    E = pow(Ec*Ec + Es*Es, 0.5);
    dt = dt - floor(En*dt/2.0/3.14159265358979) * (2.0*3.14159265358979)/En;
    Y = En*dt-Es;
    sigma = sign(Es*cos(Y) + Ec*sin(Y));
    X = Y + sigma * 0.85*E;
    s = X/pow(alpha, 0.5);
  }
  else {
    //hyperbolic motion
    A = mu/alpha;
    En = pow(-mu/A/A/A, 0.5);
    Ec = 1.0 - r0/A;
    Es = u/pow(-A*mu, 0.5);
    E = pow(Ec*Ec - Es*Es, 0.5);
    Dm = En*dt;
    if (Dm<0) s = -log((-2.0 * Dm + 1.8 * E)/(Ec - Es))/pow(-alpha, 0.5);
    else s = log((2.0 * Dm + 1.8 * E)/(Ec + Es))/pow(-alpha, 0.5);
  }
  return s;
}

int kepler_solve (mpfr_t dt, 
		  mpfr_t *F, mpfr_t *G, mpfr_t *Fdot, mpfr_t *Gdot) {
  double f, fp;
  double ds, s;
  double r0, v0s, u, alpha;
  double c0, c1, c2, c3;
  int i=0;

  r0 = pow(x*x +  y*y + z*z, 0.5);
  v0s = vx * vx + vy*vy + vz * vz;
  u = x*vx + y*vy + z*vz;
  alpha = 2.0*mu/r0 - v0s;

  s = initial_guess_for_s (dt, r0, u, alpha);
  ds = 1.0;
  while (norm(ds) > 1.0e-12) {
    stumpff(s*s*alpha, &c0, &c1, &c2, &c3);
    c1 *= s; c2 *= s*s; c3 *= s*s*s;
    f   = r0 * c1 + u * c2 + mu * c3 - dt;
    fp  = r0 * c0 + u * c1 + mu * c2;
    ds = f/fp;
    s -= ds;
    if (i++>50) {
      printf("Convergence error in Newton method\n");
      return -1; 
    }
  }
  *F = 1.0 - (mu/r0) * c2;
  *G = dt - mu * c3;
  *Fdot = - (mu/fp/r0) * c1;
  *Gdot = 1.0 - (mu/fp) * c2;
  return 0;
}

int evolve (double time_new) {
  mpfr_t x_new, y_new, z_new;
  mpfr_t vx_new, vy_new, vz_new;
  mpfr_t F, G, Fdot, Gdot;
  mpfr_t dt;
  
  mpfr_inits2(PR, 
	      x_new, y_new, z_new,
	      vx_new, vy_new, vz_new,
	      F, G, Fdot, Gdot,
	      dt,
	      (mpfr_ptr)0);

  dt = time_new - time;

  if (kepler_solve(dt, &F, &G, &Fdot, &Gdot) == 0) {  
    x_new = x * F + vx * G;
    y_new = y * F + vy * G;
    z_new = z * F + vz * G;
    vx_new = x * Fdot + vx * Gdot;
    vy_new = y * Fdot + vy * Gdot;
    vz_new = z * Fdot + vz * Gdot;
    x = x_new; y = y_new; z = z_new; vx = vx_new; vy = vy_new; vz = vz_new;
    time = time_new;
    return 0;
  }
  else return -1;

  mpfr_clears(x_new, y_new, z_new,
	      vx_new, vy_new, vz_new,
	      F, G, Fdot, Gdot,
	      dt,
	      (mpfr_ptr)0);
}

int set_mu(double mu_) {
  mu = mu_;
  return 0;
}

int set_position(double r[3]) {
  x = r[0]; y = r[1]; z = r[2];
  return 0;
}

int get_position(double r[3]) {
  r[0] = x; r[1] = y; r[2] = z;
  return 0;
}

int set_velocity(double v[3]) {
  vx = v[0]; vy = v[1]; vz = v[2];
  return 0;
}
