/*
Kepler, two body integrator.

    Copyright (C) 2010 Marcell Marosvolgyi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

-- 
    This is code demonstrates different integrators for the two-body
    problem. The methods are:
    1) Exact. 
    2) Implicit euler.
    3) Explicit euler.
    
    The differential equation being solved:
    ..
    r  = -mu * r / |r|^3
    where, m_1 r_1 + m_2 r_2 = 0 and
    r = r1 - r2
    
*/

#include <stdio.h>
#include <math.h>
#include <mpfr.h>

#define EULER_EXPLICIT 0
#define EULER_IMPLICIT 1
#define EXACT 2

double time=0.0;
double q1, q2, p1, p2;
double dt = 0.01;
double H0 = -0.5;

void initial_conditions(double e) {
  //for excentric ell. orbit..
  q1 = 1 - e;
  q2 = 0.0;
  p1 = 0.0;
  p2 = pow((1+e)/(1-e), 0.5);
  H0 = -0.5;
}

double energy() {
  //return current value of hamiltonian
  return 0.5 * (p1*p1 + p2*p2) - pow(q1*q1+q2*q2, -0.5);
}

double energy_error() {
  return H0-energy();
}

void evolve(double t_new, short mode) {
  double r3;
  double q1n, q2n;
  int steps;
  int i;

  steps = (t_new - time) / dt;

  for (i=0; i < steps; i++) {
    if (mode == EULER_IMPLICIT) {
      q1 += p1 * dt;
      q2 += p2 * dt;
      r3 = pow(q1*q1+q2*q2, 1.5); 
      p1 += -q1/r3 * dt;
      p2 += -q2/r3 * dt;
    }
    if (mode == EULER_EXPLICIT) {
      q1n = q1 + p1 * dt;
      q2n = q2 + p2 * dt;
      r3 = pow(q1*q1+q2*q2, 1.5); 
      p1 += -q1/r3 * dt;
      p2 += -q2/r3 * dt;
      q1 = q1n;
      q2 = q2n;
    }
    if (mode == EXACT) {
      
    }
  }
  time = t_new;
}

int main(int argc, char *argv[]) {
  double t;
  initial_conditions(0.7);//set excentricity i.c.
  
  for (t = 0.0; t<16.0*2*3.14159265358979; t+=0.1) {
    evolve(t, EULER_EXPLICIT);
    printf("%2.3e, %2.3e, %2.5e\n", q1, q2, energy_error());
  } 
  return 0;
}
