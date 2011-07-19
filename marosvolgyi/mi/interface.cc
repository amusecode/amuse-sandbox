#include<iostream>
#include<fstream>
#include<cmath>
#include<unistd.h>
#ifdef SAPPORO_GRAPE

#include "sapporo.h"
#include "../lib/g6/g6lib.h"
#else
//#include"../lib/g6/g6_dummy.h"
#endif //GRAPE6

//#include "src/IO.h"

//#include<mpi_interface.h>
//#include"drive.h"

#include "worker_code.h"
//#include "../lib/stopcond/stopcond.h"
#include "Particle.h"
#include "Nbody_amuse.h"

using namespace std;

extern int echo(int input);

//int echo_int(int input, int * output){
//    *output = echo(input);
//    return 0;
//}

//add and remove particles
int new_particle(int * index_of_the_particle, double mass, double radius, 
  double x, double y, double z, double vx, double vy, double vz)
{
  //need to finish
  int i;
  i = myglobals.Ntot++;
  myglobals.prt[i].mass = mass; 
  *index_of_the_particle = i;
  return 0;
}

int delete_particle(int index_of_particle)
{
  //i = index_of_particle;
  //prt[i].mass = 0.0;
  //prt[i].radius = 0.0;
  return 0;
}

//set and get index of particles
int get_index_of_first_particle(int * index_of_the_particle)
{
  //*index_of_the_particle = prt.front();
  return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle)
{
  /*
  i = index_of_the_particle;
  if (i < (prt.size() - 1))
  {
    *index_of_the_next_particle = prt.at(i+1);
    return 0;
  }
  else
  {
    if (i == (prt.size() -1))
    {
      return 1;
    }
    else
    {
      return -1;
    }
  }
  */
  return 0;
}

int get_indices_of_colliding_particles(int * index_of_particle1, 
  int * index_of_particle2)
{
  return 0;
}

//set and get: mass, total mass etc.
int set_mass(int index_of_the_particle, double setMass)
{
  /*
  if(mpi_rank == 0)     
    { // calculate only on the root mpi process, not on others
        prt[i].mass = setMass;
        return 0;
	}*/
  return 0;
}

int get_mass(int index_of_the_particle, double * getMass)
{
  /*
  i = index_of_the_particle;
  *getMass = prt[i].mass;*/
  return 0;
}

int get_total_mass(double * getMass)
{
  /*
    if(mpi_rank == 0)     
    { // calculate only on the root mpi process, not on others
      *getMass = 0.0;

      for (int i = 0; i< Ntot; i++)
      {
        *getMass += mass[i];
      }
        return 0;
	}*/
  return 0;
}

//get and set: position, velocity, and acceleration
int set_position(int index_of_the_particle, double x, double y, double z)
{/*
  i = index_of_the_particle;
  prt[i].pos[0] = x;
  prt[i].pos[1] = y;
  prt[i].pos[2] = z;*/
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
  double * z)
{/*
  i = index_of_the_particle;
  *x = prt[i].pos[0];
  *y = prt[i].pos[1];
  *z = prt[i].pos[2];*/
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
  double vz)
{/*
  i = index_of_the_particle;
  prt[i].veli[0] = vx;
  prt[i].veli[1] = vy;
  prt[i].veli[2] = vz;*/
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
  double * vz)
{/*
  i = index_of_the_particle;
  *vx = prt[i].veli[0];
  *vy = prt[i].veli[1];
  *vz = prt[i].veli[2];*/
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az)
{
  /*
  i = index_of_the_particle;
  prt[i].acci[0] = ax;
  prt[i].acci[1] = ay;
  prt[i].acci[2] = az;*/
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
  double * az)
{/*
  i = index_of_the_particle;
  *ax = prt[i].acci[0];
  *ay = prt[i].acci[1];
  *az = prt[i].acci[2];*/
  return 0;
}

//set and get: different times
int set_time(double time)
{
  return 0;
}

int get_time(double * time)
{
  return 0;
}

int set_dt_dia(double dt_dia)
{//not implimented
  return 0;
}

int get_dt_dia(double * dt_dia)
{//not implimented
  return 0;
}

int get_time_step(double * time_step)
{//not implimented
  return 0;
}

int set_dt_param(double dt_dia)
{//not implimented
  return 0;
}

int get_dt_param(double * dt_dia)
{//not implimented
  return 0;
}

//center of mass values
int get_center_of_mass_position(double * x, double * y, double * z)
{/*
  if (mpi_rank == 0)
    { // calculate only on the root mpi process, not on others
      *x=0; *y=0; *z=0;
      double m = 0.0;
      get_total_mass(&m)
  

      for (int i = 0; i<Ntot; i++)
      {
	 *x += prt[i].mass*prt[i].pos[0];
	 *y += prt[i].mass*prt[i].pos[1];
	 *z += prt[i].mass*prt[i].pos[2];
      }

      *x /= m;
      *y /= m;
      *z /= m;
      return 0;
      }*/
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{/*
  if (mpi_rank == 0)
    { // calculate only on the root mpi process, not on others
      *vx=0; *vy=0; *vz=0;
      double m = 0.0;
      get_total_mass(&m)
  

      for (int i = 0; i<Ntot; i++)
      {
         *x += prt[i].mass*prt[i].veli[0];
	 *y += prt[i].mass*prt[i].veli[1];
	 *z += prt[i].mass*prt[i].veli[2];
      }

      *vx /= m;
      *vy /= m;
      *vz /= m;
      return 0;
      }*/
  return 0;
}

//energy and gravity values
int get_potential_at_point(double eps, double x, double y, double z, 
  double * phi)
{/*
  if (mpi_rank == 0)
    { // calculate only on the root mpi process, not on others
      double r,rx,ry,rz;
      *V = 0.0;
      for (int i = 0; i < Ntot; i++)
	{
          rx = prt[i].pos[0]-x;
          ry = prt[i].pos[1]-y;
          rz = prt[i].pos[2]-z;
          r = sqrt(rx*rx+ry*ry+rz*rz + eps2);
          *phi -= prt[i].mass/r;
        }
	}*/
  return 0;
}

int get_potential(double x, double y, double z, double * V)
{
  return 0;
}

int get_kinetic_energy(double * kinetic_energy)
{
  return 0;
}

int get_potential_energy(double * potential_energy)
{
  return 0;
}

int get_gravity_at_point(double eps, double x, double y, double z, 
  double * forcex, double * forcey, double * forcez)
{/*
  if(mpi_rank == 0)   
  { // calculate only on the root mpi process, not on others
    double rx, ry, rz, r3, r2, r, F;

    *forcex = 0;*forcey = 0;*forcez = 0;

    for (int i = 0; i<Ntot; i++)
    {
        rx = prt[i].pos[0] - x;
        ry = prt[i].pos[1] - y;
        rz = prt[i].pos[2] - z;
        r2 = (rx*rx+ry*ry+rz*rz + eps2);//do I need to make cases for each of eps2_fs_smbh, eps2_fs_imbh, eps2_bh,
        r = sqrt(r2);
        r3 = r2*r;
        F = prt[i].mass/r3;
        *forcex += F * rx;
        *forcey += F * ry;
        *forcez += F * rz;
    }
    return 0;
    }*/
  return 0;
}

//misc
int set_pair_detect_factor(double pair_detect_factor)
{
  return 0;
}

int get_pair_detect_factor(double * pair_detect_factor)
{
  return 0;
}

int get_number_of_particles(int * number_of_particles)
{
  //*number_of_particles = Ntot;
  return 0;
}

int get_radius(int index_of_the_particle, double * Thisradius)
{
  //i = index_of_the_particle;
  //*Thisradius = prt[i].radius;
  return 0;
}

int get_total_radius(double * radius)
{
  // not implemented
  return 0;
}

int set_eps2(double epsilon_squared)//need to think about these two...do I need seperat cases for bh vs fs particles?
{
  printf("Don't use set_eps2 to set this parameter as, for this code, eps2 varies in this code depending on the type of interaction (ie. if it is star-star interaction or star-BH interaction).  Use instead set_smoothing.");
  return 0;
}

int get_eps2(double * epsilon_squared)
{
  printf("Don't use get_eps2 to get this parameter as, for this code, eps2 varies in this code depending on the type of interaction (ie. if it is star-star interaction or star-BH interaction).  Use instead get_smoothing.");
  return 0;
}

//set and get: state of a given particle
int set_state(int index_of_the_particle, double mass, double radius, 
  double x, double y, double z, double vx, double vy, double vz)
{/*
  if (mpi_rank == 0)
  { // calculate only on the root mpi process, not on others
    i = index_of_the_particle;
    if (i > prt.size())
      {
        return -1;
      }
    mass = prt[i].mass;
    radius = prt[i].radius;
    prt[i].pos[0] = x;
    prt[i].pos[1] = y;
    prt[i].pos[2] = z;
    prt[i].veli[0] = vx;
    prt[i].veli[1] = vy;
    prt[i].veli[2] = vz;
    return 0;
    }*/
  return 0;
}

int get_state(int index_of_the_particle, double * mass, double * radius, 
  double * x, double * y, double * z, double * vx, double * vy, 
  double * vz)
{/*
  if (mpi_rank == 0)
  { // calculate only on the root mpi process, not on others
    i = index_of_the_particle;
    if (i > prt.size())
      {
        return -1;
      }
    *mass = prt[i].mass;
    *radius = prt[i].radius;
    *x = prt[i].pos[0];
    *y = prt[i].pos[1];
    *z = prt[i].pos[2];
    *vx = prt[i].veli[0];
    *vy = prt[i].veli[1];
    *vz = prt[i].veli[2];
    return 0;
    }*/
  return 0;
}


//get various stopping condition information
int get_stopping_condition_info(int index, int * type, 
  int * number_of_particles)
{
  return 0;
}

int get_stopping_condition_particle_index(int index, 
  int index_of_the_column, int * index_of_particle)
{
  return 0;
}

int set_stopping_condition_out_of_box_parameter(double value)
{
  return 0;
}

int get_stopping_condition_out_of_box_parameter(double * value)
{
  return 0;
}

int get_stopping_condition_number_of_steps_parameter(int * value)
{
  return 0;
}

int set_stopping_condition_number_of_steps_parameter(int value)
{
  return 0;
}

int get_stopping_condition_timeout_parameter(double * value)
{
  return 0;
}

int get_number_of_stopping_conditions_set(int * result)
{
  return 0;
}

int set_stopping_condition_timeout_parameter(double value)
{
  return 0;
}

//non-set/get interactions
int evolve(double time)
{
  return 0;
}

int evolve_model(double time)
{
  return 0;
}

int get_potential(int id, double *phi)
{
  return -1;
}

int set_radius(int id, double r)
{
  return -1;
}

int get_smoothing(int id1, int id2, double *eps)
{
  return 0;
}

int set_smoothing(int id1, int id2, double eps)
{
  return 0;
}


int initialize_code()
{
  return 0;
}

int synchronize_model()
{
  return 0;
}

int enable_stopping_condition(int type)
{
  return 0;
}

int is_stopping_condition_set(int type, int * result)
{
  return 0;
}

int is_stopping_condition_enabled(int type, int * result)
{
  return 0;
}

int disable_stopping_condition(int type)
{
  return 0;
}

int has_stopping_condition(int type, int * result)
{
  return 0;
}

int cleanup_code()
{
  return 0;
}

int recommit_particles()
{//do I need these?
  return 0;
}

int commit_particles()
{
  return 0;
}

int recommit_parameters()
{
  return 0;
}

int commit_parameters()
{
  return 0;
}

