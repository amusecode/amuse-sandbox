#include "cvode/cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode/cvode_dense.h"				/* prototype for CVDense */
#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "cvode/sundials_types.h"			/* definition of type realtype */

int evolve(
    double m1, double m2, double m3,
    double R1, double R2, double R3,
    double a1, double a2,
    double e1, double e2,
    double INCL1, double INCL2, double AP1, double AP2, double LAN1, double LAN2,
    double t, double dt,
    double * m1_out, double * m2_out, double * m3_out,
    double * R1_out, double * R2_out, double * R3_out,
    double * a1_out, double * a2_out,
    double * e1_out, double * e2_out,
    double *INCL1_out, double *INCL2_out, double *INCL12_out, double * AP1_out, double * AP2_out, double *LAN1_out, double *LAN2_out,
    double * t_out,
    int * out_flag, int * error_flag);
    
static int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data);
static int fev_triad(realtype t, N_Vector yev, N_Vector ydot, void *data);
static int froot(realtype t, N_Vector yev, realtype *gout, void *data);	/* ODE Rootfinding function */

int get_equations_of_motion_specification(int *relative_tolerance_t);
int set_equations_of_motion_specification(int relative_tolerance_t);
int get_relative_tolerance(double *relative_tolerance_t);
int set_relative_tolerance(double relative_tolerance_t);
int get_f_oct(double *f_oct_t);
int set_f_oct(double f_oct_t);
int get_f_1PN_1(double *f_1PN_1_t);
int set_f_1PN_1(double f_1PN_1_t);
int get_f_1PN_2(double *f_1PN_2_t);
int set_f_1PN_2(double f_1PN_2_t);
int get_f_1PN_12(double *f_1PN_12_t);
int set_f_1PN_12(double f_1PN_12_t);
int get_f_25PN_1(double *f_25PN_1_t);
int set_f_25PN_1(double f_25PN_1_t);
int get_f_25PN_2(double *f_25PN_2_t);
int set_f_25PN_2(double f_25PN_2_t);

int compute_triad_vectors_from_orbital_elements(double m1, double m2, double m3, double a1, double a2, double e1, double e2, double INCL1, double INCL2, double AP1, double AP2, double LAN1, double LAN2,
    double e1_vec[3], double e2_vec[3], double h1_vec[3], double h2_vec[3], double q1_vec[3], double q2_vec[3]);
int compute_orbital_elements_from_triad_vectors(double m1, double m2, double m3,
    double e1_vec[3],double e2_vec[3],double h1_vec[3],double h2_vec[3],double q1_vec[3],double q2_vec[3],
    double *a1_out, double *a2_out, double *e1_out, double *e2_out, double *INCL1_out, double *INCL2_out, double *AP1_out, double *AP2_out, double *LAN1_out, double *LAN2_out);
    
void cross3(double a[3], double b[3], double result[3]);
double norm3(double v[3]);
double norm3_squared(double v[3]);
double dot3(double a[3], double b[3]);
