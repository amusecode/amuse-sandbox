//prototypes
int get_center_of_mass_velocity(double *cmvx, double *cmvy, double *cmvz);
int plot(int R, int G, int B);
int initialization(int precision);
int get_precision(int *precision);
int viewer(int view_);
int set_position(int id_, double x_, double y_, double z_);
int set_velocity(int id_, double vx_, double vy_, double vz_);
int get_position(int id_, double *x_, double *y_, double *z_);
int get_velocity(int id_, double *vx_, double *vy_, double *vz_);
int get_velocity_mp(int id_, char *xmp, char *ymp, char *zmp);
int set_reduced_mass(double mu_);
int commit_particles();
int set_mass(int id_, double mass);
int get_center_of_mass(double *cmx, double *cmy, double *cmz);
int get_center_of_mass_velocity(double *cmvx, double *cmvy, double *cmvz);
int evolve_system(double new_time);
int get_number_of_particles(int *number_of_particles);
int new_particle(int *id_, double mass_, double radius, 
		 double x_, double y_, double z_,
                 double vx_, double vy_, double vz_);

double midX = Xres/2;
double midY = Yres/2;

struct {
  double x0, y0, z0;
  double x1, y1, z1;
  double vx0, vy0, vz0;
  double vx1, vy1, vz1;
  //char vx0mp[500];
  //char vy0mp[500];
  //char vz0mp[500];
  //char vx1mp[500];
  //char vy1mp[500];
  //char vz1mp[500];
  mpfr_t vx0mp;
  mpfr_t vy0mp;
  mpfr_t vz0mp;
  mpfr_t vx1mp;
  mpfr_t vy1mp;
  mpfr_t vz1mp;
  double mass1;
  double mass2;
  double radius1;
  double radius2;
#ifdef __PLOT
  SDL_Surface *screen;
#endif
  short view;
} my_globals;
