#include "cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode_dense.h"				/* prototype for CVDense */
#include "sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "sundials_types.h"			/* definition of type realtype */

int evolve(double m1, double m2, double m3, double R1, double R2, double a1, double a2, double e1, double e2, double itot, double g1, double g2, double t, double dt, double * m1_out, double * m2_out, double * m3_out, double * R1_out, double * R2_out, double * a1_out, double * a2_out, double * e1_out, double * e2_out, double * itot_out, double * g1_out, double * g2_out, double * t_out, int * out_flag, int * error_flag);
static int fev(realtype t, N_Vector yev, N_Vector ydot, void *data);
static int froot(realtype t, N_Vector yev, realtype *gout, void *data);	/* ODE Rootfinding function */
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

