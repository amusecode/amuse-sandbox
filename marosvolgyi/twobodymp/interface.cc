/* 
   AMUSE interface to kepler.c
   ===========================
   
   This code also interfaces two particles to c.o.m. system
   of kepler.c. For description see Fundamentals of Celestial 
   Mechanics, J.M.A. Danby 2nd Edition, chapter 6 introduction.

 */

#define __PLOT
#include "worker_code.h"
#include "src/kepler.h"

#ifdef __PLOT
#ifdef __cplusplus
extern "C" {
#endif
#include <SDL/SDL.h>
#include "../include/draw.h"
#ifdef __cplusplus
}
#endif

#define Xres 1000
#define Yres 800

#endif
#include "interface.h"

int plot(int R, int G, int B) {
#ifdef __PLOT
  double x_= 0.0;
  double y_= 0.0;
  double z_= 0.0;
  int rR, rG, rB;
  int mR, mG, mB;

  get_position(0, &x_, &y_, &z_);
  DrawCircle(my_globals.screen, 200.0*x_+midX, 200.0*y_+midY, my_globals.radius1, 10, R, G, B);
  if (id == 1) {
    get_position(1, &x_, &y_, &z_);
    DrawCircle(my_globals.screen, 200.0*x_+midX, 200.0*y_+midY, my_globals.radius2, 10, R, G, B);
  }
  SDL_Flip(my_globals.screen);
#endif
  return 0;
}

int viewer(int view_) {
  if ((my_globals.view == 0) & (view_ == 1)) {
    my_globals.view = view_;
#ifdef __PLOT
    SDL_Init( SDL_INIT_EVERYTHING );
    my_globals.screen = SDL_SetVideoMode( Xres, Yres, 32, SDL_SWSURFACE );
    my_globals.view = 1;
    DrawText(my_globals.screen, "KEPLER", 1, 1, 's', 0,0,255);
#endif
  }
  return 0;
}

int initialization(int precision) {
  initialize(precision);
  id = -1;
  my_globals.view = 0;
  my_globals.mass1 = 0.0;
  my_globals.mass2 = 0.0;

  return 0;
}

int get_precision(int *precision) {
  *precision = PR;
  return 0;
}

int set_position(int id_, double x_, double y_, double z_) {
  if (id_ == 0) {
    my_globals.x0 = x_; my_globals.y0 = y_; my_globals.z0 = z_;
  }
  if (id_ == 1) {
    my_globals.x1 = x_; my_globals.y1 = y_; my_globals.z1 = z_;
  }

  if (id == 0) {
    set_pos(x_, y_, z_);
    return 0;
  }
  else if (id == 1) {
    double dx = my_globals.x0 - my_globals.x1;
    double dy = my_globals.y0 - my_globals.y1;   
    double dz = my_globals.z0 - my_globals.z1;
    set_pos(dx, dy, dz);
    return 0;
  }
  return -1;
}

int get_position(int id_, double *x_, double *y_, double *z_) {
  if (id == 0) get_pos(x_, y_, z_);
  if (id == 1) {
    double factor;
    if (id_ == 0) factor = (my_globals.mass1/my_globals.mass2 + 1);
    else if (id_ == 1) factor = -(my_globals.mass2/my_globals.mass1 + 1);
    get_pos(x_, y_, z_);      
    *x_ /= factor;
    *y_ /= factor;
    *z_ /= factor;
    double cmx =0.0;
    double cmy =0.0;
    double cmz =0.0;

    get_center_of_mass(&cmx, &cmy, &cmz);
    *x_ += cmx;
    *y_ += cmy;
    *z_ += cmz;
  }    
  return 0;
}

int set_velocity(int id_, double vx_, double vy_, double vz_) {
  if (id_ == 0 ) {
    my_globals.vx0 = vx_; my_globals.vy0 = vy_; my_globals.vz0 = vz_;
  }
  if (id_ == 1) {
    my_globals.vx1 = vx_; my_globals.vy1 = vy_; my_globals.vz1 = vz_;
  }
  
  if (id == 0) {
    set_vel(vx_, vy_, vz_);
    return 0;
  }
  else if (id ==1 ) {
    double dvx = my_globals.vx0 - my_globals.vx1;
    double dvy = my_globals.vy0 - my_globals.vy1;   
    double dvz = my_globals.vz0 - my_globals.vz1;
    set_vel(dvx, dvy, dvz);
    return 0;
  }
  return -1;
}

int get_velocity(int id_, double *vx_, double *vy_, double *vz_) {
  if (id == 0) get_vel(vx_, vy_, vz_);
  if (id == 1) {
    double factor;
    if (id_ == 0) factor = (my_globals.mass1/my_globals.mass2 + 1);
    else if (id_ == 1) factor = -(my_globals.mass2/my_globals.mass1 + 1);
    get_vel(vx_, vy_, vz_);      
    *vx_ /= factor;
    *vy_ /= factor;
    *vz_ /= factor;
    double cmvx =0.0;
    double cmvy =0.0;
    double cmvz =0.0;

    get_center_of_mass_velocity(&cmvx, &cmvy, &cmvz);
    *vx_ += cmvx;
    *vy_ += cmvy;
    *vz_ += cmvz;
  }    
  return 0;
}

int set_reduced_mass(double mu_) {
  set_mu(mu_);
  return 0;
}
int set_mass(int id_, double mass) {
  if (id_ ==  0) my_globals.mass1 = mass;
  if (id_ ==  1) my_globals.mass2 = mass;
  return 0;
}


int set_radius(int id_, double radius) {
  if (id_ == 0) my_globals.radius1 = radius;
  if (id_ == 1) my_globals.radius2 = radius;
  return 0;
}
   
int get_center_of_mass(double *cmx, double *cmy, double *cmz) {
  if (id == 1 ) {
    *cmx = my_globals.x0 * my_globals.mass1 +
      my_globals.x1 * my_globals.mass2; 
    *cmy = my_globals.y0 * my_globals.mass1 +
      my_globals.y1 * my_globals.mass2; 
    *cmz = my_globals.z0 * my_globals.mass1 +
      my_globals.z1 * my_globals.mass2; 
  }
  return 0;
}

int get_center_of_mass_velocity(double *cmvx, double *cmvy, double *cmvz) {
  if (id == 1 ) {
    *cmvx = my_globals.vx0 * my_globals.mass1 +
      my_globals.vx1 * my_globals.mass2; 
    *cmvy = my_globals.vy0 * my_globals.mass1 +
      my_globals.vy1 * my_globals.mass2; 
    *cmvz = my_globals.vz0 * my_globals.mass1 +
      my_globals.vz1 * my_globals.mass2; 
  }
  return 0;
}

int get_kinetic_energy(int id_, double *Ek) {
  double vx_;
  double vy_;
  double vz_;
  double mass_;

  if (id == 1) {
    get_velocity(id_, &vx_, &vy_, &vz_);
    if (id_ == 0) mass_ = my_globals.mass1;
    if (id_ == 1) mass_ = my_globals.mass2;
    *Ek = 0.5 * mass_ * vx_ * vx_ + vy_ * vy_ + vz_ * vz_;

    return 0;
  }
  return -1;
}

int new_particle(int *id_, double mass_, double radius, 
		 double x_, double y_, double z_,
                 double vx_, double vy_, double vz_) {
  
  if (id == -1) {
    // define single particle description
    *id_ = 0;
    id = 0;
  }    
  else if (id == 0) {
    *id_ = 1;
    id = 1;
  }
  set_position(id, x_, y_, z_);
  set_velocity(id, vx_, vy_, vz_);
  set_mass(id, mass_);
  set_radius(id, radius);
  return 0;
}

int get_number_of_particles(int *number_of_particles) {
  *number_of_particles = id;
  return 0;
}

int commit_particles() {
  //asume two are set
  if (id == 1) {
    set_reduced_mass(my_globals.mass1 + my_globals.mass2);
    return 0;
  }
  else if (id == 0) {
    set_reduced_mass(my_globals.mass1);
  }
  return -1;

}

int evolve_system(double new_time) {
  evolve_d(new_time);
  return 0;
}

