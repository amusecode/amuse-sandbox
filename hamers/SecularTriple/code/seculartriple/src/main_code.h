#include <math.h>
#include <cstdlib>

#include "cvode/cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode/cvode_dense.h"				/* prototype for CVDense */
#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "cvode/sundials_types.h"			/* definition of type realtype */

#include "helper_routines.h"
#include "ODE_system.h"

#define CONST_G			    (double)	6.67259e-8
#define CONST_G_P2          (double)    CONST_G*CONST_G
#define CONST_G_P3          (double)    CONST_G*CONST_G_P2
#define CONST_C_LIGHT		(double)	2.99792458e10
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define c_1div2             (double)    1.0/2.0
#define c_3div2             (double)    3.0/2.0
#define c_5div2             (double)    5.0/2.0
#define c_15div2            (double)    15.0/2.0
#define c_31div2            (double)    31.0/2.0
#define c_1div4             (double)    1.0/4.0
#define c_15div4            (double)    15.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_3div5             (double)    3.0/5.0
#define c_8div5             (double)    8.0/5.0
#define c_64div5            (double)    64.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_3div8             (double)    3.0/8.0
#define c_15div8            (double)    15.0/8.0
#define c_45div8            (double)    45.0/8.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0
#define c_1div16            (double)    1.0/16.0
#define c_5div16            (double)    5.0/16.0
#define c_9div16            (double)    9.0/16.0
#define c_15div16           (double)    15.0/16.0
#define c_185div16          (double)    185.0/16.0
#define c_11div18           (double)    11.0/18.0
#define c_73div24           (double)    73.0/24.0
#define c_5div64            (double)    5.0/64.0
#define c_25div64           (double)    25.0/64.0
#define c_37div96           (double)    37.0/96.0
#define c_121div304         (double)    121.0/304.0

/*	ODE solver related quantities	*/
#define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef struct {
	double m1,m2,m3;
    double R1,R2,R3;
    double m1dot;
    double AMC_star1,AMC_star2,AMC_star3;
    double gyration_radius_star1,gyration_radius_star2,gyration_radius_star3;
    double k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3;
    double gamma_in,gamma_out,mu_in,mu_out;
    int f_quad,f_oct,f_tides,f_mass_transfer,f_1PN_in,f_1PN_out,f_1PN_in_out,f_25PN_in,f_25PN_out;
} *UserData;

int evolve(
    double m1, double m2, double m3,
    double R1, double R2, double R3,
    double spin_angular_frequency1, double spin_angular_frequency2, double spin_angular_frequency3,
    double AMC_star1, double AMC_star2, double AMC_star3,
    double gyration_radius_star1, double gyration_radius_star2, double gyration_radius_star3,
    double k_div_T_tides_star1, double k_div_T_tides_star2, double k_div_T_tides_star3,
    double a_in, double a_out,
    double e_in, double e_out,
    double INCL_in, double INCL_out, double AP_in, double AP_out, double LAN_in, double LAN_out,
    double m1dot,
    double gamma_in, double gamma_out, double mu_in, double mu_out,       
    double t, double dt,
    double * m1_output, double * m2_output, double * m3_output,
    double * R1_output, double * R2_output, double * R3_output,
    double * spin_angular_frequency1_output, double * spin_angular_frequency2_output, double * spin_angular_frequency3_output,
    double * a_in_output, double * a_out_output,
    double * e_in_output, double * e_out_output,
    double *INCL_in_output, double *INCL_out_output, double *INCL_in_out_output, double * AP_in_output, double * AP_out_output, double *LAN_in_output, double *LAN_out_output,
    double * t_output,
    int * output_flag, int * error_flag);
    
int get_equations_of_motion_specification(int *relative_tolerance_t);
int set_equations_of_motion_specification(int relative_tolerance_t);
int get_relative_tolerance(double *relative_tolerance_t);
int set_relative_tolerance(double relative_tolerance_t);
int get_f_quad(double *f_quad_t);
int set_f_quad(double f_quad_t);
int get_f_oct(double *f_oct_t);
int set_f_oct(double f_oct_t);
int get_f_mass_transfer(double *f_mass_transfer_t);
int set_f_mass_transfer(double f_mass_transfer_t);
int get_f_tides(double *f_tides_t);
int set_f_tides(double f_tides_t);
int get_f_1PN_in(double *f_1PN_in_t);
int set_f_1PN_in(double f_1PN_in_t);
int get_f_1PN_out(double *f_1PN_out_t);
int set_f_1PN_out(double f_1PN_out_t);
int get_f_1PN_in_out(double *f_1PN_in_out_t);
int set_f_1PN_in_out(double f_1PN_in_out_t);
int get_f_25PN_in(double *f_25PN_in_t);
int set_f_25PN_in(double f_25PN_in_t);
int get_f_25PN_out(double *f_25PN_out_t);
int set_f_25PN_out(double f_25PN_out_t);
