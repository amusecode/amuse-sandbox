/* 
Kepler integrator, two-body problem

code based on:
Fundamentals of Celestial Mechanics, J.M.A. Danby 2nd Edition
*/

#include <stdio.h>
#include <math.h>

#define DEBUG 0
#define DEBUG_S 1

double x, y, z;
double vx, vy, vz;
double time=0;
double mu;

int factorial (int k) {
  int r = 1;
  int i;

  for (i=1; i<=k; i++) {
    r *=i;
  }
  return r;
}

double norm(double x) {
  if (x<0) return -x;
  return x;
}

double sign(double x) {
  if (norm(x)<1e-12) return 0.0;
  if (x<0) return -1.0;
  else return 1.0;
}

void stumpff(double s, double *c0, double *c1, double *c2, double *c3) {
  double sqrt_s;

  if (s>0.0) {
    sqrt_s = pow(s, 0.5);
    *c0 = cos(sqrt_s);
    *c1 = sin(sqrt_s)/sqrt_s;
  }
  else if (s<0.0){
    sqrt_s = pow(-s, 0.5);
    *c0 = cosh(sqrt_s);
    *c1 = sinh(sqrt_s)/sqrt_s;
  }
  else printf("Error in stumpff s = 0\n");
  *c2 = 1.0/s * (1.0 - *c0);
  *c3 = 1.0/s * (1.0 - *c1);
}

double initial_guess_for_s (double dt, 
			    double r0, double u, double alpha,
			    double *c0, double *c1, double *c2, double *c3) {
  double A, En, Ec, Es, E, X, Y, Dm, sigma;
  double s;

  if (norm(dt/r0) <=.2) {
    if (DEBUG) printf("other..\n");
    s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0);
  }
  else if (alpha>0) {
    //elliptic motion initial guess
    if (DEBUG) printf("elliptic\n");
    A = mu/alpha;
    En = pow(mu/A/A/A, 0.5);
    Ec = 1.0 - r0/A;
    Es = u/En/A/A;
    E = pow(Ec*Ec + Es*Es, 0.5);
    dt = dt - floor(En*dt/2/3.14159265358979) * (2*3.14159265358979)/En;
    Y = En*dt-Es;
    sigma = sign(Es*cos(Y) + Ec*sin(Y));
    X = Y + sigma * 0.85*E;
    s = X/pow(alpha, 0.5);
  }
  else {
    //hyperbolic motion
    if (DEBUG) printf("hyperbolic\n");
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

int keplers(double dt, 
	     double r0, double u, double alpha, 
	     double *fp_, 
	     double *c0, double *c1, double *c2, double *c3) {
  double f, fp, fpp, fppp;
  double ds, s;
  int i=0;

  s = initial_guess_for_s (dt, r0, u, alpha, c0, c1, c2, c3);
  if (DEBUG) printf("s_guess=%2.3e\n", s);
  ds = 1.0;
  while (norm(ds)>1.0e-12) {
    stumpff(s*s*alpha, c0, c1, c2, c3);
    if (DEBUG) printf("c = %2.3e,c = %2.3e, c=%2.3e\n", *c0, *c1, *c2);
    *c1 *= s; *c2 *= s*s; *c3 *= s*s*s;
    f   = r0 * *c1 + u * *c2 + mu * *c3 - dt;
    fp  = r0 * *c0 + u * *c1 + mu * *c2;
    //*fpp = (-r0*alpha + mu) * c1 - u * c0;
    //*fppp= (-r0*alpha + mu) * c0 - u * alpha * c1;
    ds = f/fp;
    s -= ds;
    if (DEBUG) printf("f = %2.3e, fp = %2.3e, ds=%2.3e\n", f, fp, ds);
    if (i++>150) {
      printf("Convergence error in Newton method\n");
      return -1; 
    }
  }
  if (DEBUG) printf("s_final=%2.3e, fp=%2.3e\n", s, fp);
  *fp_ = fp;
  return 0;
}

int evolve(double time_new) {
  double r0, v0s, u, alpha;
  double fp, c0, c1, c2, c3;
  double dt;
  double F, G, Fdot, Gdot;
  double x_new, y_new, z_new;
  double vx_new, vy_new, vz_new;

  r0 = pow(x*x +  y*y + z*z, 0.5);
  v0s = vx * vx + vy*vy + vz * vz;
  u = x*vx + y*vy + z*vz;
  alpha = 2.0*mu/r0 - v0s;
  dt = time_new - time;

  if (keplers(dt, r0, u, alpha, &fp, &c0, &c1, &c2, &c3) == 0) {  
    F = 1.0 - (mu/r0)* c2;
    G = dt - mu * c3;
    Fdot = - (mu/fp/r0) * c1;
    Gdot = 1.0 - (mu/fp) * c2;
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
}

int main(int argc, char *argv[]) {
  double t;
  FILE *datafile;
  mu = 1.0;
  
  datafile = fopen("kepler.dat","w");
  x = 1.0; y = 0.1; z = -0.1; vx = -0.1; vy = 2.0; vz = -0.2;
  for (t = 0.01; t<5.0; t+=0.05) {
    if (evolve(t)==-1) {
      //printf("Newton foot finding failed\n");
      //break;
    }
    fprintf(datafile, "%2.3e, %2.3e, %2.3e\n", x, y, z);
  }
  fflush(datafile);
  fclose(datafile);
  return 0;
}

/*

*/
