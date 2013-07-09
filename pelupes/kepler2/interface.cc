#include "src/stdinc.h"
#include "src/kepler.h"

#include <map>

static double simtime;
static map<int, kepler*> kmap;
static int kcounter=0;

static double central_mass=0.;

using namespace std;

int initialize_code()
{
    return 0;
}

int commit_parameters()
{
    return 0;
}

int commit_particles()
{
    return 0;
}

int recommit_particles()
{
    return 0;
}

int synchronize_model()
{
    return 0;
}

int recommit_parameters()
{
    return 0;
}

int cleanup_code()
{
    return 0;
}

int set_acceleration(int particle_identifier, double ax, double ay, double az) {
    return -2; // Not implemented
}

int get_acceleration(int particle_identifier, double *ax, double *ay, double *az) {
    return -2; // Not implemented
}
int get_potential(int id,double *p){
    return -2; // Not implemented
}


int set_central_mass(double m){ central_mass=m;return 0;}
int get_central_mass(double *m){ *m=central_mass;return 0;}

int set_eps2(double m){ return 0;}
int get_eps2(double *m){ *m=0.;return 0;}

int get_time_step(double *m){ *m=0.;return 0;}


int new_particle(int *id, double mass,
			double x, double y, double z,
			double vx, double vy, double vz,double r)
{
    kepler *k = new kepler;
    k->set_time(simtime);
    k->set_total_mass(mass + central_mass);
    k->set_rel_pos(vec(x,y,z));
    k->set_rel_vel(vec(vx,vy,vz));
    k->initialize_from_pos_and_vel();
    *id=kcounter++;
   	kmap.insert(kmap.end(), std::pair<int, kepler*>( *id, k));

    return 0;
}

int delete_particle(int id)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    
    if (iter != kmap.end()){
        delete (*iter).second;
        kmap.erase(iter);
    } else {
      return -3;
    }
    return 0;
}

int get_state(int id, double *m, double *x, double *y, double *z,
                 double *vx, double *vy, double *vz,double *rad)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    
    if (iter != kmap.end()){
        kepler *k=iter->second;
        *m=k->get_total_mass()-central_mass;
        vec r = k->get_rel_pos();
        *x = r[0];
        *y = r[1];
        *z = r[2];
        vec v = k->get_rel_vel();
        *vx = v[0];
        *vy = v[1];
        *vz = v[2];
        *rad = 0.;
    } else {
      return -3;
    }
    return 0;
}

int set_state(int id, double m, double x, double y, double z,
                 double vx, double vy, double vz,double r)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        kepler *k=iter->second;
        k->set_total_mass(m+central_mass);
        k->set_rel_pos(vec(x,y,z));
        k->set_rel_vel(vec(vx,vy,vz));
        k->initialize_from_pos_and_vel();
    } else {
      return -3;
    }
    return 0;
}


int set_mass(int id, double m)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        kepler *k=iter->second;
        k->set_total_mass(m+central_mass);
        k->initialize_from_pos_and_vel();
    } else {
      return -3;
    }
    return 0;
}

int set_radius(int id, double r)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
    } else {
      return -3;
    }
    return 0;
}

int get_radius(int id, double *r)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
      *r=0;
    } else {
      return -3;
    }
    return 0;
}

int set_position(int id, double x, double y, double z)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        kepler *k=iter->second;
        k->set_rel_pos(vec(x,y,z));
        k->initialize_from_pos_and_vel();
    } else {
      return -3;
    }
    return 0;
}

int set_velocity(int id,double vx, double vy, double vz)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        kepler *k=iter->second;
        k->set_rel_vel(vec(vx,vy,vz));
        k->initialize_from_pos_and_vel();
    } else {
      return -3;
    }
    return 0;
}

int evolve_model(double tend)
{
   map<int, kepler*>::iterator iter;
//   for(iter=kmap.begin(); iter!=kmap.end(); iter++)
//   {
//      iter->second->transform_to_time(tend);       
//   }
   int keys[kmap.size()];
   int i=0;
   for(iter=kmap.begin(); iter!=kmap.end(); iter++)
   {
      keys[i]=iter->first;
      i++;
   }
#pragma omp parallel for   
   for(i=0;i<kmap.size();i++)
   {
      kmap[keys[i]]->transform_to_time(tend);       
   }

   simtime=tend;
   return 0;
}

int get_mass(int id, double *mass)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       *mass=k->get_total_mass()-central_mass;
    } else {
      return -3;
    }
    return 0;
}

int get_time(double * t)
{
    *t = simtime;
    return 0;
}

int get_position(int id,double * x, double * y, double * z)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       vec r = k->get_rel_pos();
       *x = r[0];
       *y = r[1];
       *z = r[2];
    } else {
      return -3;
    }
    return 0;
}

int get_velocity(int id, double * vx, double * vy, double * vz)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       vec v = k->get_rel_vel();
       *vx = v[0];
       *vy = v[1];
       *vz = v[2];
    } else {
      return -3;
    }
    return 0;
}

int get_semi_major_axis(int id,double * semi)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       *semi = k->get_semi_major_axis();
    } else {
      return -3;
    }
    return 0;
}

int get_eccentricity(int id,double * ecc)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       *ecc = k->get_eccentricity();
    } else {
      return -3;
    }
    return 0;
}

int get_specific_orbital_energy(int id,double * ecc)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       *ecc = k->get_energy();
    } else {
      return -3;
    }
    return 0;
}


int get_next_radial_crossing_time(int id, double radius, double *tnext)
{
    map<int, kepler*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       kepler *k=iter->second;
       if( k->get_energy() < 0) 
       {
         if(radius > k->get_apastron() || radius < k->get_periastron() )return -5;
         *tnext=k->pred_advance_to_radius(radius);
       } else 
       {
         if(radius < k->get_periastron() ) return -6;
         if(radius > k->get_separation() && k->get_true_anomaly() > 0) return -7;          
         *tnext=k->pred_advance_to_radius(radius);
         if(*tnext < k->get_time()) return -8;
       }
    } else {
      return -3;
    }
    return 0;

 
}



int get_index_of_first_particle(int *particle_identifier) {
    return -2; // Not implemented
}
int get_index_of_next_particle(int particle_identifier, int *next_particle_identifier) {
    return -2; // Not implemented
}

int get_indices_of_colliding_particles(int *index_of_particle1, int *index_of_particle2) {
    return -2; // Not implemented
}

// simulation property getters:
int get_total_mass(double *total_mass) {
    return -2; // Not implemented
}
int get_total_radius(double *total_radius) {
    return -2; // Not implemented
}

int set_begin_time(double input) {
    return -2; // Not implemented
}

int get_begin_time(double * output) {
    return -2; // Not implemented
}

int get_center_of_mass_position(double *x, double *y, double *z){
    return -2; // Not implemented
}
int get_center_of_mass_velocity(double *vx, double *vy, double *vz) {
    return -2; // Not implemented
}
int get_kinetic_energy(double *kinetic_energy) {
    return -2; // Not implemented
}
int get_potential_energy(double *potential_energy) {
    return -2; // Not implemented
}
int get_number_of_particles(int *number_of_particles) {
    *number_of_particles = kmap.size();
    return 0;
}


