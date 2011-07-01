#define FLOAT  double
#define CLFLOAT cl_double
#define CLFLOAT4 cl_double4
#define DOUBLE long double
#define FLT(x) __CONCAT(x,)
#define DBL(x) __CONCAT(x,l)
#define INT int
#define UINT unsigned int
#define LONG long
#define ULONG unsigned long

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}

#define LOG(fmt, ...) {\
	printf("%s:%d\t", __FILE__, __LINE__);\
	printf(fmt, ## __VA_ARGS__);\
}

#define ABS(X) (((X) >= 0) ? (X) : -(X))

#define ENDRUN(fmt, ...) { \
	printf("ENDRUN@%s:%d ", __FILE__, __LINE__);\
	printf(fmt, ## __VA_ARGS__);\
	fflush(stdout);\
	exit(-1);\
}

struct particle {
  UINT id;
  FLOAT mass;
  FLOAT radius; /*not used*/
  DOUBLE pos[3];
  DOUBLE vel[3];
  DOUBLE pot;
  DOUBLE postime;
  FLOAT timestep;
  INT level; 

  DOUBLE kick_amount;
  //UINT visited_cc;
};

struct jparticle {
  FLOAT mass;
  FLOAT pos[3];
  FLOAT vel[3];
};

struct sys {
  UINT n; 
  struct particle *part;
  struct particle *last;
  struct sys *next_cc;
};

enum intopt {
  CONSTANT, //0
  UNSPLIT,  //1
  PASS,     //2
  HOLD,     //3
  BRIDGE,   //4
  NAIVE,    //5
  VARIABLE, //6
  PASS_DKD, //7
  HOLD_DKD, //8
  CC_SPLIT2, //9
  UNSPLIT4, //10
  CC_SPLIT4,  //11
  EVOLVE_OK2, //12
  EVOLVE_OK4, //13
  TWOBODY, //14
  CC_SPLIT2_TWOBODY, //15
  EVOLVE_ROK2 //16
};

extern FLOAT eps2;
extern FLOAT dt_param;

void init_code();
void stop_code();
void init_evolve(struct sys s);
void do_evolve(struct sys s, double dt, int inttype);
DOUBLE system_potential_energy(struct sys s);
DOUBLE system_kinetic_energy(struct sys s);

//void get_evolve_statistics_(int *ttot, int *ktot, int *dtot, int *tstot, int *kstot, int *dstot);
void get_evolve_statistics_(double *ttot, double *ktot, double *dtot, double *tstot, double *kstot, double *dstot, double *cetot, double *cetotfail);
DOUBLE sys_initial_timestep(struct sys s);
DOUBLE sys_forces_min_timestep(struct sys s);
DOUBLE ok_timestep_ij_fw(struct particle *i, struct particle *j);

DOUBLE system_linear_momentum_x(struct sys s);
DOUBLE system_center_of_mass_x(struct sys s);
DOUBLE system_angular_momentum_x(struct sys s);
