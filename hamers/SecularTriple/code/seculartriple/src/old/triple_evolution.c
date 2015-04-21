#include "../binary_prototypes.h"
#include "../binary_structures.h"
#include "../binary_macros.h"
#include "../iteration/iteration_macros.h"
#include <math.h>

#include "triple_prototypes.h"
#include "../maths_triple/cvode.h"				/* prototypes for CVODE fcts., consts. */
#include "../maths_triple/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "../maths_triple/cvode_dense.h"			/* prototype for CVDense */
#include "../maths_triple/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "../maths_triple/sundials_types.h"			/* definition of type realtype */

/*	Frequently used constants; GRAVITATIONAL_CONSTANT and M_SUN are defined in binary_macros.h	*/
#define CONST_G			(double)	GRAVITATIONAL_CONSTANT
#define CONST_GM		(double)	CONST_G*M_SUN
#define CONST_GMM		(double)	CONST_GM*M_SUN
#define CONST_GM_SQRT		(double)	sqrt(CONST_GM)
#define CONST_GM_POW3		(double)	pow(CONST_GM,3.0)
#define CONST_ANG_TRIPLE	(double)	M_SUN*CONST_GM_SQRT
#define CONST_C_LIGHT		(double)	2.99792458e10	/* For now not defined as SPEED_OF_LIGHT because of the numerical error (see binary_macros.h) */		
#define CONST_C_LIGHT_POW2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_DIV_POW_5	(double)	pow(CONST_C_LIGHT,-5.0)	
#define CONST_GR		(double)	CONST_GM_POW3*CONST_C_LIGHT_DIV_POW_5
	
#define TIDES_CRIT_STAR1 ( (fabs(e1dot_third_component) < fabs(e1dot_tides_star1)) && (fabs(e1dot_third_component) > 0.1*fabs(e1dot_tides_star1)) )
#define TIDES_CRIT_STAR2 ( (fabs(e1dot_third_component) < fabs(e1dot_tides_star2)) && (fabs(e1dot_third_component) > 0.1*fabs(e1dot_tides_star2)) )
#define GWE_CRIT ( (fabs(e1dot_third_component) < fabs(e1dot_GR_diss)) && (fabs(e1dot_third_component) > 0.1*fabs(e1dot_GR_diss)) )


/*	ODE solver related quantities	*/
#define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define NEQ		8				/* Number of ODE equations */
#define RTOL		RCONST(1.0e-10)			/* Scalar relative tolerance - acceptable value is 1.0e-10 as determined by trial and error */
#define ATOL1		RCONST(1.0e-8)
#define ATOL2		RCONST(1.0e-8)
#define ATOL3		RCONST(1.0e-8)
#define ATOL4		RCONST(1.0e-8)
#define ATOL5		RCONST(1.0e1)
#define ATOL6		RCONST(1.0e-8)
#define ATOL7		RCONST(1.0e-8)
#define ATOL8		RCONST(1.0e-8)				
#define TIMESTEP0	RCONST(10.0*YEAR_LENGTH_IN_SECONDS)		/* Initial internal ODE timestep */
#define MAXNUMSTEPS 	5e8						/* Maximum number of internal steps */
#define MAXTIME		RCONST(13.7e9*YEAR_LENGTH_IN_SECONDS)		/* Maximum integration time	*/
#define MAXNUMCONVFAIL	20						/* Maximum number of convergence failures in fev */

//#define A1_CHANGE_TOL		(double)	0.1
//#define OMSPIN1_CHANGE_TOL	(double)	1.0
//#define OMSPIN2_CHANGE_TOL	(double)	1.0

//#define R_CHANGE_TOL		(double)	0.00001
//#define M_CHANGE_TOL		(double)	0.00000001
//#define KT_CHANGE_TOL		(double)	0.00001
//#define L_CHANGE_TOL		(double)	0.001

/*	The following structure is used to pass on user data from the main function to the ODE right hand side function fev	*/
typedef struct {
	struct stardata_t * stardata_u;
	double a1_initial,a2_initial;						/* Unit: cm */
	double m1,m2,m3;							/* Unit: M_SUN */
	double radius_star1,radius_star2;					/* Unit: cm */
	double apsidal_motion_constant_star1,apsidal_motion_constant_star2;
	double gyration_radius_star1,gyration_radius_star2;
	double kT_star1, kT_star2;
	double omega_spin_star1_initial,omega_spin_star2_initial;
	double t_start;
	double e1_initial;
	int ST1,ST2;
	int thread_id;
	Boolean RLOF_boolean;
} *UserData;

static int fev(realtype t, N_Vector yev, N_Vector ydot, void *data);			/* ODE right-hand-side */
static int froot(realtype t, N_Vector yev, realtype *gout, void *data);			/* ODE Rootfinding function */
static void ehfun(int error_code, const char *module, const char *function, char *msg, void *eh_data);
static void PrintFinalStats(void *cvode_mem);						/* Function for ODE solver printout */
static int check_flag(void *flagvalue, char *funcname, int opt);			/* Function used by ODE solver */

void triple_evolution(struct stardata_t * RESTRICT stardata, Boolean RLOF_boolean)
{

	if ((stardata->star[1].radius*pow(stardata->star[1].rol,-1.0) > 0.99) || (stardata->star[2].radius*pow(stardata->star[2].rol,-1.0) > 0.99)) return;

	struct common_t *common = &(stardata->common);	
	struct star_t * star1 = &(stardata->star[1]);
	struct star_t * star2 = &(stardata->star[2]);

//	common->triple_log_double_array[0] = stardata->model.model_time;
//	common->triple_log_double_array[1] = stardata->common.separation/AURSUN;
//	common->triple_log_double_array[2] = stardata->common.triple_outer_semi_major_axis;
//	common->triple_log_double_array[3] = log10(1.0-stardata->common.eccentricity);
//	common->triple_log_double_array[4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);
//	common->triple_log_double_array[5] = common->triple_total_inclination*180.0/M_PI;
//	common->triple_log_double_array[6] = stardata->star[1].mass;
//	common->triple_log_double_array[7] = stardata->star[2].mass;
//	common->triple_log_double_array[8] = stardata->common.triple_outer_mass;
//	common->triple_log_double_array[9] = stardata->star[1].radius;
//	common->triple_log_double_array[10] = stardata->star[2].radius;
//	common->triple_log_int_array[0] = 0;
//	common->triple_log_int_array[1] = stardata->star[1].stellar_type;
//	common->triple_log_int_array[2] = stardata->star[2].stellar_type;

/*	Note: the "common->triple_log_double_array" variables are now assigned in iteration/main_loop.c	*/

//	if ((common->triple_active == TRUE) && (stardata->model.model_time > 215)) return;

	Boolean change_timestep = FALSE;
//	Boolean output_log;
//	output_log = TRUE;			/* Set this to TRUE to enable triple updates in the output file */

//	common->triple_new_timestep = 0.0;	/* If zero, then the timestep does not have to be modified on the account of triple_c */

	/*	The initial values of e1, e2, g1, g2 and the total inclination - angles in units of rad	*/
	double e1_initial = common->eccentricity;
	double e2_initial = common->triple_outer_orbit_eccentricity;
	double g1_initial = common->triple_inner_orbit_argument_of_periastron;
	double g2_initial = common->triple_outer_orbit_argument_of_periastron;
	double total_inclination_initial = common->triple_total_inclination;
	if (total_inclination_initial == 0) total_inclination_initial = 0.01*M_PI/180.0; /* If cos(itot) = 1, then the solver will have trouble with the triple ODEs, so instead set the mutual inclination angle to a small value */
//	if (total_inclination_initial == 180) total_inclination_initial = 180.1;

	/*	Spin angular frequencies in binary_c are in units of rad/year; in triple_c the units are rad/s	*/
	double omega_spin_star1_initial = stardata->star[1].omega/YEAR_LENGTH_IN_SECONDS;	/* Unit: rad/s */
	double omega_spin_star2_initial = stardata->star[2].omega/YEAR_LENGTH_IN_SECONDS;	/* Unit: rad/s */

	double m1 = stardata->star[1].mass;			/* Unit: M_SUN */
	double m2 = stardata->star[2].mass;			/* Unit: M_SUN */
	double m3 = common->triple_outer_mass;			/* Unit: M_SUN */

	double radius_star1,radius_star2,I_star1,I_star2;
	if (RLOF_boolean == FALSE)
	{
		I_star1 = moment_of_inertia(star1,stardata->star[1].radius);	/* Unit: Msun Rsun^2 */
		I_star2 = moment_of_inertia(star2,stardata->star[2].radius);	/* Unit: Msun Rsun^2 */
		radius_star1 = R_SUN*stardata->star[1].radius;	/* Units: cm */
		radius_star2 = R_SUN*stardata->star[2].radius;	/* Units: cm */
	}
	else
	{
		I_star1 = moment_of_inertia(star1,stardata->star[1].radx);	/* Unit: Msun Rsun^2 */
		I_star2 = moment_of_inertia(star2,stardata->star[2].radx);	/* Unit: Msun Rsun^2 */
		radius_star1 = R_SUN*stardata->star[1].radx;	/* Units: cm */
		radius_star2 = R_SUN*stardata->star[2].radx;	/* Units: cm */
	}

	/*	Frequently used mass combinations (for speedup)	*/
	double m_prod_inner = m1*m2;
	double m_tot_inner = m1+m2;
	double m_tot_triple = m_tot_inner+m3;

	/*	The initial values of a1 and a2; note that by this point changes to a2 due to mass loss in the inner/outer binary systems have already been taken into account by the triple_c function triple_calc_mass_transfer_effects_on_outer_orbit */
	double P1_initial = common->orbital_period; 			/* Unit: years */
	double P2_initial = common->triple_outer_orbital_period; 	/* Unit: years */

	double a1_initial = R_SUN*common->separation;				/* Unit: cm */
	double a2_initial = R_SUN*AURSUN*common->triple_outer_semi_major_axis;	/* Unit: cm */

	/*	Compute critical ratio (a2/a1)_c for which the triple system is hierarchically stable	(Marding Aarseth (2001), Eq. 90))	*/
	double beta_initial = a2_initial/a1_initial;
	double q_out = m3/m_tot_inner;
	double beta_crit = (1.0/(1.0-e2_initial))*2.8*pow( (1.0+q_out)*(1.0+e2_initial)/sqrt(1.0-e2_initial),2.0/5.0)*(1.0 \
		- 0.3*total_inclination_initial/M_PI);
	if (beta_initial < beta_crit) 
	{
		if (common->triple_cml_logging == TRUE)	
		{
			if (common->triple_log_cases[4] == 0)
			{
				printf("The current ratio a2/a1 is %e, which is smaller than the critical value of %e. \n Hence the triple system is no longer dynamically stable! Skipping triple evolution. \n",beta_initial,beta_crit);

				FILE *fp_cml;
				char buff_cml[50] = "";
				sprintf(buff_cml,"triple_cml_output%d.txt",stardata->common.triple_thread_id);
				fp_cml = fopen(buff_cml,"a");
				/* Output variables: CASE t/Myr SNT a1/AU a2/AU e1 e2 itot/deg ST1 ST2 M1 M2 R1 R2 */
				fprintf(fp_cml,"4\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\t%g\t%g\t%g\t%g\n",stardata->model.model_time,0,stardata->common.separation/AURSUN,stardata->common.triple_outer_semi_major_axis,log10(1.0-stardata->common.eccentricity),log10(1.0-stardata->common.triple_outer_orbit_eccentricity),stardata->common.triple_total_inclination*180.0/M_PI,stardata->star[1].stellar_type,stardata->star[2].stellar_type,stardata->star[1].mass,stardata->star[2].mass,stardata->common.triple_outer_mass,stardata->star[1].radius,stardata->star[2].radius);
				fclose(fp_cml);
				common->triple_log_cases[4] += 1;
			}
		}
		return;
	}

	/*	Compute typical Kozai cycle timescale	*/
	double Pkozai = 4.0*pow(beta_initial,3.0/2.0)*sqrt( ( pow(a1_initial*beta_initial,3.0)*m_tot_triple ) / (CONST_GM) )*(1.0/m3)*pow(1-e2_initial*e2_initial,3.0/2.0);

//	if ((stardata->model.model_time > 150) && (stardata->model.model_time < 250))
//	{
//		stardata->model.dtm = 0.1;
//	}

	/*	Compute the initial value of Gatot: total angular momentum of the triple system	*/
	double L1_const = CONST_ANG_TRIPLE*m_prod_inner*sqrt(1.0/m_tot_inner);		/* Constant terms in L1 */
	double L1_initial = L1_const*sqrt(a1_initial);					/* Units: g cm^2 s^-1 */	
	double L2_const = CONST_ANG_TRIPLE*m3*m_tot_inner*sqrt(1.0/m_tot_triple);	/* Constant terms in L2 */
	double L2_initial = L2_const*sqrt(a2_initial);					/* Units: g cm^2 s^-1 */	
	double Ga1_initial = L1_initial*sqrt(1.0 - e1_initial*e1_initial);		
	double Ga2_initial = L2_initial*sqrt(1.0 - e2_initial*e2_initial);		
	double Gatot_initial = sqrt(pow(Ga1_initial,2.0) + pow(Ga2_initial,2.0) \
		+ 2.0*Ga1_initial*Ga2_initial*cos(total_inclination_initial));		

	/*	Obtain current time and binary_c timestep - all in units of seconds	*/
	double t_start = stardata->model.model_time*1.0e6*YEAR_LENGTH_IN_SECONDS;	/* Note: stardata->model.model_time: unit: Myr */
	double timestep_binary_c = stardata->model.dtm*1.0e6*YEAR_LENGTH_IN_SECONDS;	/* Note: stardata->model.dtm: unit: Myr */
	double t_end = t_start + timestep_binary_c;

	double t_reached,timestep_new,cvode_last_timestep;
	t_reached = timestep_new = cvode_last_timestep = 0.0;

	/*	Some stellar parameters required for tides, assumed to be constant during triple integration	*/
//	double apsidal_motion_constant_star1 = 0.1*ENVELOPE_MASS(1)/((double) m1);	/* From Hurley's thesis, p. 80 */
//	double apsidal_motion_constant_star2 = 0.1*ENVELOPE_MASS(2)/((double) m2);	/* From Hurley's thesis, p. 80 */
	double apsidal_motion_constant_star1 = common->triple_amc_star1;	/* Calculated previously in triple_calc_apsidal_motion_constants */
	double apsidal_motion_constant_star2 = common->triple_amc_star2;	/* Calculated previously in triple_calc_apsidal_motion_constants */
	double gyration_radius_star1,gyration_radius_star2;
	if (RLOF_boolean == FALSE)
	{
		gyration_radius_star1 = I_star1/(stardata->star[1].mass*stardata->star[1].radius*stardata->star[1].radius);
		gyration_radius_star2 = I_star2/(stardata->star[2].mass*stardata->star[2].radius*stardata->star[2].radius);
	}
	else
	{
		gyration_radius_star1 = I_star1/(stardata->star[1].mass*stardata->star[1].radx*stardata->star[1].radx);
		gyration_radius_star2 = I_star2/(stardata->star[2].mass*stardata->star[2].radx*stardata->star[2].radx);
	}
	
	double factor_kT1 = 1.0;
	double factor_kT2 = 1.0;

	/* AD HOC KT MULTIPLICATION */
//	factor_kT1 = 1000.0;
//	factor_kT2 = 1000.0;
//	if (stardata->star[1].stellar_type <= 1) {factor_kT1 = 1000.0;}
//	if (stardata->star[2].stellar_type <= 1) {factor_kT2 = 1000.0;}

	double kT_star1 = factor_kT1*triple_tidal_constants(1,stardata);	/* The constant (k/T) in Hut's equations, for star 1 */
	double kT_star2 = factor_kT2*triple_tidal_constants(2,stardata);	/* The constant (k/T) in Hut's equations, for star 2 */

//	if (stardata->star[1].radius*pow(stardata->star[1].rol,-1.0) > 0.99)
//	{
//		return 0;
//	}
	
//	if (stardata->star[2].radius*pow(stardata->star[2].rol,-1.0) > 0.99)
//	{
//		return 0;
//	}
	
	/*	ODE right hand side function constants	*/
	UserData data;
	UserData eh_data;
	data = NULL;
	data = (UserData) malloc(sizeof *data);
	data->stardata_u = stardata;
	data->a1_initial = a1_initial;
	data->a2_initial = a2_initial;
	data->m1 = m1;
	data->m2 = m2;
	data->m3 = m3;
	data->radius_star1 = radius_star1;
	data->radius_star2 = radius_star2;
	data->apsidal_motion_constant_star1 = apsidal_motion_constant_star1;
	data->apsidal_motion_constant_star2 = apsidal_motion_constant_star2;
	data->gyration_radius_star1 = gyration_radius_star1;
	data->gyration_radius_star2 = gyration_radius_star2;
	data->kT_star1 = kT_star1;
	data->kT_star2 = kT_star2;
	data->omega_spin_star1_initial = omega_spin_star1_initial;
	data->omega_spin_star2_initial = omega_spin_star2_initial;
	data->t_start = t_start;
	data->e1_initial = e1_initial;
	data->RLOF_boolean = RLOF_boolean;
	data->ST1 = stardata->star[1].stellar_type;
	data->ST2 = stardata->star[2].stellar_type;
	data->thread_id = stardata->common.triple_thread_id;

	eh_data = data;

	/*	Initialize ODE solver variables	*/
	N_Vector yev, yev_reached, abstol;
	void *cvode_mem;
	int flag,flag_s;
	int nroot = 2;					/* Number of roots to be found by rootfinding routine */
	int rootsfound[nroot];				/* Integer array which will contain information on which roots were found */

	yev = yev_reached = abstol = NULL;
	cvode_mem = NULL;
	
	yev = N_VNew_Serial(NEQ);
	if (check_flag((void *)yev, "N_VNew_Serial", 0)) return;
	yev_reached = N_VNew_Serial(NEQ);
	if (check_flag((void *)yev_reached, "N_VNew_Serial", 0)) return;
	abstol = N_VNew_Serial(NEQ); 
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return;
	
		
	/*	Initialize y-vector; to avoid singularities in the evolution equations, e1 and e2 need to be set to very small values	*/
	if (e1_initial <= 0.0) {
		e1_initial = 1e-10;
	}
	if (e2_initial <= 0.0) {
		e2_initial = 1e-10;
	}

//	Ith(yev,1) = e1_initial;		/* Inner orbit eccentricity */
//	Ith(yev,2) = e2_initial;		/* Outer orbit eccentricity */	
	Ith(yev,1) = log10(1.0 - e1_initial);
	Ith(yev,2) = log10(1.0 - e2_initial);
	Ith(yev,3) = g1_initial;		/* Inner orbit argument of periastron - unit: rad */
	Ith(yev,4) = g2_initial;		/* Outer orbit argument of periastron - unit: rad */
	Ith(yev,5) = a1_initial;		/* Inner orbit semi-major axis - unit: cm */
//	Ith(yev,6) = Gatot_initial;		/* Total angular momentum - unit: g cm^2 s^-1 (contains information on inclination angle) */
	Ith(yev,6) = cos(total_inclination_initial);
	Ith(yev,7) = omega_spin_star1_initial;	/* Spin angular frequency of star 1 - unit: rad s^-1 */
	Ith(yev,8) = omega_spin_star2_initial;	/* Spin angular frequency of star 2 - unit: rad s^-1 */
	
	/*	Set relative and absolute tolerances	*/
	Ith(abstol,1) = ATOL1;   
	Ith(abstol,2) = ATOL2;
	Ith(abstol,3) = ATOL3;
	Ith(abstol,4) = ATOL4;
	Ith(abstol,5) = ATOL1;   
	Ith(abstol,6) = ATOL2;
	Ith(abstol,7) = ATOL3;
	Ith(abstol,8) = ATOL4;
	
	/*	Initialize CVode	*/
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return;
	
	flag = CVodeInit(cvode_mem, fev, 0.0, yev);
//	flag = CVodeInit(cvode_mem, fev, t_start, yev);
	if (check_flag(&flag, "CVodeInit", 1)) return;

	flag = CVodeSetErrHandlerFn(cvode_mem, ehfun, eh_data);
	if (check_flag(&flag, "CVodeSetErrHandlerFn", 1)) return;

	flag = CVodeSVtolerances(cvode_mem, RTOL, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return;
	
	flag = CVDense(cvode_mem, NEQ);
	if (check_flag(&flag, "CVDense", 1)) return;
  
	flag = CVodeSetInitStep(cvode_mem, TIMESTEP0);
	if (check_flag(&flag, "CVodeSetInitStep", 1)) return;
		
	flag = CVodeSetMaxNumSteps(cvode_mem, MAXNUMSTEPS);
	if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return;

//	flag = CVodeSetMinStep(cvode_mem, 0.1);
//	if (check_flag(&flag, "CVodeSetMinStep", 1)) return;

	flag = CVodeSetMaxHnilWarns(cvode_mem, 1);
	if (check_flag(&flag, "CVodeSetMaxHnilWarns", 1)) return;
			
//	flag = CVodeSetStopTime(cvode_mem, MAXTIME);
//	if (check_flag(&flag, "CVodeSetStopTime", 1)) return(1);

	flag = CVodeSetMaxConvFails(cvode_mem, MAXNUMCONVFAIL);
	if (check_flag(&flag, "CVodeSetMaxConvFails", 1)) return;
			
	flag = CVodeSetUserData(cvode_mem, data);
	if (check_flag(&flag, "CVodeSetUsetData", 1)) return;

	flag = CVodeRootInit(cvode_mem, nroot, froot);
	if (check_flag(&flag, "CVodeRootInit", 1)) return;	

	/* ============================================================================================================
	* Run the solver - new y-vector is yev_reached
	*/


	double t_end_cvode = t_end - t_start;
	double t_reached_cvode = 0.0;

//	flag_s = CVode(cvode_mem, t_end, yev_reached, &t_reached, CV_NORMAL);
	flag_s = CVode(cvode_mem, t_end_cvode, yev_reached, &t_reached_cvode, CV_NORMAL);	

	t_reached = t_reached_cvode + t_start;

	CVodeGetLastStep(cvode_mem, &cvode_last_timestep);

	if (flag_s == CV_ROOT_RETURN)
	{
		CVodeGetRootInfo(cvode_mem,rootsfound);
		timestep_new = t_reached - t_start;
		change_timestep = TRUE;
		stardata->model.dt = timestep_new/YEAR_LENGTH_IN_SECONDS;
		stardata->model.dtm = stardata->model.dt*1.0E-06;
		stardata->model.the_timestep = stardata->model.dt;	
	}

//	if ((flag_s != CV_SUCCESS) && (flag_s != CV_ROOT_RETURN))
//	{
//		printf("CVODE_ERROR %d id %d time %g\n",flag_s,common->triple_thread_id,stardata->model.model_time);
//		FILE *tempfile;
//		tempfile = fopen("CVODE_errors.txt","a");
//		fprintf(tempfile,"error %d id %d time %g\n",flag_s,common->triple_thread_id,stardata->model.model_time);
//		fclose(tempfile);
//	}
		
//	if (rootsfound[0] == 1) printf("COLLISION during triple integration \n");
//	if (rootsfound[0] == -1) printf("COLLISION during triple integration \n");
//	if (rootsfound[1] == 1) printf("ang mom change dt; new dt = %g \n",timestep_new*1.0e-6/YEAR_LENGTH_IN_SECONDS);
//	if (rootsfound[1] == -1) printf("ang mom change dt; new dt = %g \n",timestep_new*1.0e-6/YEAR_LENGTH_IN_SECONDS);

	/* ============================================================================================================
	* The new values
	*/

//	double e1_reached = Ith(yev_reached,1);
//	double e2_reached = Ith(yev_reached,2);

	double x_reached = Ith(yev_reached,1);
	double y_reached = Ith(yev_reached,2);	
	double e1_reached = 1.0 - pow(10,x_reached);
	double e2_reached = 1.0 - pow(10,y_reached);
	double g1_reached = Ith(yev_reached,3);
	double g2_reached = Ith(yev_reached,4);
	double a1_reached = Ith(yev_reached,5);
//	double Gatot_reached = Ith(yev_reached,6);
	double cositot_reached = Ith(yev_reached,6);
	if (cositot_reached > 1.0)
	{
		cositot_reached = 2.0 - cositot_reached;
	}

	if (cositot_reached < -1.0)
	{
		cositot_reached = -2.0 - cositot_reached;
	}

	double total_inclination_reached = acos(cositot_reached);
	double omega_spin_star1_reached = Ith(yev_reached,7);
	double omega_spin_star2_reached = Ith(yev_reached,8);
	
	double a2_reached = a2_initial;				/* Assume a2 changes only due to isotropic wind in the inner orbit - this is taken care of by the function triple_calc_mass_transfer_effects_on_outer_orbit */
	double L1_reached = L1_const*sqrt(a1_reached);
	double L2_reached = L2_const*sqrt(a2_reached);
	double Ga1_reached = L1_reached*sqrt(1.0 - e1_reached*e1_reached);
	double Ga2_reached = L2_reached*sqrt(1.0 - e2_reached*e2_reached);
#ifdef IGNORE
	double cos_total_inclination_reached = (pow(Gatot_reached,2.0) - pow(Ga1_reached,2.0) \
		- pow(Ga2_reached,2.0))/(2.0*Ga1_reached*Ga2_reached);
	if (cos_total_inclination_reached < -1.0) 
	{
		printf("t = %e Myr; cositot < -1.0\n",t_start/(1e6*YEAR_LENGTH_IN_SECONDS));
		cos_total_inclination_reached = -1.0;
	}
	if (cos_total_inclination_reached > 1.0) 
	{
		printf("t = %e Myr; cositot > 1.0\n",t_start/(1e6*YEAR_LENGTH_IN_SECONDS));
		cos_total_inclination_reached = 1.0;
	}

	double total_inclination_reached = acos(cos_total_inclination_reached);
#endif
	double Gatot_reached = sqrt(pow(Ga1_reached,2.0) + pow(Ga2_reached,2.0) \
		+ 2.0*Ga1_reached*Ga2_reached*cos(total_inclination_reached));	
	
	/* ============================================================================================================
	* Update the relevant values in binary_c	
	*/

	stardata->common.eccentricity = e1_reached;
	stardata->common.ecc1 = e1_reached;

	common->orbital_angular_frequency = TWOPI*sqrt((m1+m2)*pow(a1_reached/(R_SUN*AURSUN),-3.0));
	common->orbital_angular_momentum = m1*m2/(m1+m2)*sqrt(1.0-e1_reached*e1_reached)*pow(a1_reached/R_SUN,2.0)*common->orbital_angular_frequency;

//	common->orbital_angular_frequency *= sqrt(pow(a1_initial/a1_reached,3.0));
//	common->orbital_angular_momentum *= sqrt(a1_reached*(1.0 - e1_reached*e1_reached)/(a1_initial*(1.0 - e1_initial*e1_initial)));
	
	double dspint1 = (omega_spin_star1_reached - omega_spin_star1_initial)*YEAR_LENGTH_IN_SECONDS;
	double dspint2 = (omega_spin_star2_reached - omega_spin_star2_initial)*YEAR_LENGTH_IN_SECONDS;


	stardata->star[1].omega += dspint1;
	stardata->star[2].omega += dspint2;

	stardata->star[1].jspin = I_star1*stardata->star[1].omega;
	stardata->star[2].jspin = I_star2*stardata->star[2].omega;


	/*	These quantities are intrinsic to triple_c	*/
	common->triple_outer_orbit_eccentricity = e2_reached;
	common->triple_inner_orbit_argument_of_periastron = g1_reached;
	common->triple_outer_orbit_argument_of_periastron = g2_reached;
	common->triple_total_inclination = total_inclination_reached;
	
	common->triple_outer_semi_major_axis = a2_reached/(R_SUN*AURSUN);
//	common->triple_outer_orbital_period = P2_initial*pow(a2_reached/a2_initial_bms,3.0/2.0)*pow((m3+m_inner_old)/(m3+m_tot_inner),1.0/2.0);
	common->triple_outer_orbital_period = pow(pow(common->triple_outer_semi_major_axis,3.0)/(m_tot_inner+m3),1.0/2.0);	/* Unit: yr */


	/*	Print some values if debug_mode is true	*/
	if ((common->triple_debug_mode == TRUE) && (stardata->model.model_time >= common->triple_debug_start_time))
	{
		printf("============================================================================================== \n");
		printf("Triple evolution debug mode \n");
		printf("============================================================================================== \n");
		if (check_flag(&flag_s, "CVode", 1)) printf("CVode error\n");
 		printf("Triple evolution start time: %g Myr\n",t_start/(1e6*YEAR_LENGTH_IN_SECONDS));
		if (flag_s == CV_SUCCESS) 
		{
			printf("Successfully evolved to time: %g Myr \n",1e-6*t_end/YEAR_LENGTH_IN_SECONDS);
			printf("Timestep: %g Myr\n",timestep_binary_c/(1e6*YEAR_LENGTH_IN_SECONDS)); 
		}
		if (flag_s == CV_TSTOP_RETURN) printf("Maximum time reached. \n");
		if (flag_s == CV_ROOT_RETURN)
		{
			if ((rootsfound[0] == 1) || (rootsfound[0] == -1)) printf("COLLISSION at t = %g Myr during triple integration - triple_c will change timestep such that the next time is precisely the time of collision. \n",t_reached/(1e6*YEAR_LENGTH_IN_SECONDS));
			if ((rootsfound[1] == 1) || (rootsfound[1] == -1)) printf("Triple_c will change the timestep from %g Myr to %g Myr such that G1 does not change more than 2 per cent\n",timestep_binary_c*1e-6/YEAR_LENGTH_IN_SECONDS,timestep_new*1e-6/YEAR_LENGTH_IN_SECONDS);
//			printf("Triple_evolution has changed binary_c timestep because during integration, the quantity ");
//			if (rootsfound[0] == -1) printf("a1 has decreased more than %g per cent. ",A1_CHANGE_TOL);
//			if (rootsfound[1] == 1) printf("a1 has increased more than %g per cent. ",A1_CHANGE_TOL);
//			if (rootsfound[2] == -1) printf("omega_spin_star1 has decreased more than ~%g per cent. ",OMSPIN1_CHANGE_TOL);
//			if (rootsfound[3] == 1) printf("omega_spin_star1 has increased more than ~%g per cent. ",OMSPIN1_CHANGE_TOL);
//			if (rootsfound[4] == -1) printf("omega_spin_star2 has decreased more than ~%g per cent. ",OMSPIN2_CHANGE_TOL);
//			if (rootsfound[5] == 1) printf("omega_spin_star2 has increased more than ~%g per cent. ",OMSPIN2_CHANGE_TOL);
//			printf("\nThe new imposed binary_c timestep is %g Myr \n",1e-6*timestep_new/YEAR_LENGTH_IN_SECONDS);
		}
//		if (common->triple_new_timestep != 0)
//		{
//			printf("Such that Ga1 does not change by more than 2 per cent, triple_c will change the timestep from %g to %g for the next iteration\n",stardata->model.dt,common->triple_new_timestep);
//		}
		printf("--------------------------------------------------------------------------------------------- \n");
		printf("M1 = %g Msun; M2 = %g Msun; M3 = %g Msun \n",m1,m2,m3);
		printf("Stellar types: k1 = %d; k2 = %d \n",stardata->star[1].stellar_type, stardata->star[2].stellar_type);
		printf("R1 = %g Rsun; RL1 = %g Rsun; R1/RL1 = %g \n",stardata->star[1].radius, stardata->star[1].rol, \ 
			stardata->star[1].radius*pow( stardata->star[1].rol,-1.0));
		printf("R2 = %g Rsun; RL2 = %g Rsun; R2/RL2 = %g \n",stardata->star[2].radius, stardata->star[2].rol, \ 
			stardata->star[2].radius*pow( stardata->star[2].rol,-1.0));
		printf("dmloss1 = %g; dmacc1 = %g \n",stardata->star[1].dmr, stardata->star[1].dmt);
		printf("dmloss2 = %g; dmacc2 = %g \n",stardata->star[2].dmr, stardata->star[2].dmt);	
		printf("k/T_1 = %g /s; k/T_2 = %g /s \n",kT_star1,kT_star2);
		printf("Gyration radii: rg1 = %g; rg2 = %g \n",gyration_radius_star1,gyration_radius_star2);
		printf("Apsidal motion constants: amc_1 = %g; amc_2 = %g \n",apsidal_motion_constant_star1,apsidal_motion_constant_star2);
		printf("Istar1 = %g; Istar2 = %g\n",I_star1,I_star2);
		printf("--------------------------------------------------------------------------------------------- \n");
		printf("Pkozai = %g Myr \n",Pkozai/(1e6*YEAR_LENGTH_IN_SECONDS));
		printf("Estimated no. of cycles: %g \n",(t_reached-t_start)/Pkozai);
		printf("--------------------------------------------------------------------------------------------- \n");
		printf("Initial values at t = %g Myr: \n",t_start/(1e6*YEAR_LENGTH_IN_SECONDS));
		printf(" P1 = %g d; P2 = %g yr; a1 = %g AU; a2 = %g AU; beta = %g; a1_periastron = %g Rsun; a1_periastron/(R1+R2) = %g\n",P1_initial*365.25,P2_initial,a1_initial/(R_SUN*AURSUN),a2_initial/(R_SUN*AURSUN),a2_initial/a1_initial,a1_initial/(R_SUN)*(1.0-e1_initial),a1_initial/R_SUN*(1.0-e1_initial)/(stardata->star[1].radius+stardata->star[2].radius));
		printf(" e1 = %g; e2 = %g; g1 = %g deg; g2 = %g deg; \n a1 = %12.10g Rsun; Gatot = %g g cm^2 s^-1; omega_spin_star1 = %g rad/s; omega_spin_star2 = %g rad/s; \n itot = %g deg \n",e1_initial,e2_initial,g1_initial*180.0/PI,g2_initial*180.0/PI, a1_initial/R_SUN, Gatot_initial, omega_spin_star1_initial,omega_spin_star2_initial,total_inclination_initial*180.0/PI);
		printf("--------------------------------------------------------------------------------------------- \n");
		printf("Final values: at t = %g Myr: \n",t_reached/(1e6*YEAR_LENGTH_IN_SECONDS));
		printf(" P1 = %g d; P2 = %g yr; a1 = %g AU; a2 = %g AU; beta = %g; a1_periastron = %g Rsun; a1_periastron/(R1+R2) = %g\n",sqrt(pow(a1_reached/(R_SUN*AURSUN),3.0)/(m1+m2))*365.25,common->triple_outer_orbital_period,a1_reached/(R_SUN*AURSUN),a2_reached/(R_SUN*AURSUN),a2_reached/a1_reached,a1_reached/(R_SUN)*(1.0-e1_reached),a1_reached/R_SUN*(1.0-e1_reached)/(stardata->star[1].radius+stardata->star[2].radius));
		printf(" e1 = %g; e2 = %g; g1 = %g deg; g2 = %g deg; \n a1 = %g Rsun; Gatot = %g g cm^2 s^-1; omega_spin_star1 = %g rad/s; omega_spin_star2 = %g rad/s; \n itot = %g deg \n",e1_reached,e2_reached,g1_reached*180.0/PI,g2_reached*180.0/PI, a1_reached/R_SUN, Gatot_reached, omega_spin_star1_reached,omega_spin_star2_reached,total_inclination_reached*180.0/PI);
		printf("--------------------------------------------------------------------------------------------- \n");	
		printf("Relative changes (in per cent): \n");
		printf(" e1: %g; e2: %g; g1: %g; g2: %g; \n a1: %g; Gatot: %g; omega_spin_star1: %g; omega_spin_star2: %g; \n itot: %g \n",100.0*(e1_reached/e1_initial-1.0),100.0*(e2_reached/e2_initial-1.0),100.0*(g1_reached/g1_initial-1.0),100.0*(g2_reached/g2_initial-1.0),100.0*(a1_reached/a1_initial-1.0),100.0*(Gatot_reached/Gatot_initial-1.0),100.0*(omega_spin_star1_reached/omega_spin_star1_initial-1.0),100.0*(omega_spin_star2_reached/omega_spin_star2_initial-1.0),100.0*(total_inclination_reached/total_inclination_initial-1.0));
		printf("--------------------------------------------------------------------------------------------- \n");	

		double e1dot_third_component_p = stardata->common.triple_e1dot_third_component;
		double e1dot_tides_star1_p = stardata->common.triple_e1dot_tides_star1;
		double e1dot_tides_star2_p = stardata->common.triple_e1dot_tides_star2;
		double e1dot_GR_diss_p = stardata->common.triple_e1dot_GR_diss;

		double a1dot_tides_star1_p = stardata->common.triple_a1dot_tides_star1;
		double a1dot_tides_star2_p = stardata->common.triple_a1dot_tides_star2;
		double a1dot_GR_diss_p = stardata->common.triple_a1dot_GR_diss;

//		double e1_p2 = e1_initial*e1_initial;
//		double e1_p4 = e1_p2*e1_p2;
//		double e1_p6 = e1_p2*e1_p4;
//		double e1_p8 = e1_p4*e1_p4;
//		double f1_p = 1.0 + (31.0/2.0)*e1_p2 + (255.0/8.0)*e1_p4 + (185.0/16.0)*e1_p6 + (25.0/64.0)*e1_p8;
//		double f2_p = 1.0 + (15.0/2.0)*e1_p2 + (45.0/8.0)*e1_p4 + (5.0/16.0)*e1_p6;
//		double f_GR_adot_p = 1.0 + (73.0/24.0)*e1_p2 + (37.0/96.0)*e1_p4;
//		double e1_p2com = 1.0 - e1_p2;
//		double omegaquot1_p = stardata->star[1].omega/stardata->common.orbital_angular_frequency;
//		double omegaquot2_p = stardata->star[2].omega/stardata->common.orbital_angular_frequency;
//		double a1dot_tides_star1_p = -6.0*(1.0+m2/m1)*kT_star1*(m2/m1)*pow(stardata->star[1].radius/stardata->common.separation,8.0)*stardata->common.separation*R_SUN*pow(e1_p2com,-15.0/2.0)*(f1_p - pow(e1_p2com,3.0/2.0)*f2_p*omegaquot1_p);
//		double a1dot_tides_star2_p = -6.0*(1.0+m1/m2)*kT_star2*(m1/m2)*pow(stardata->star[2].radius/stardata->common.separation,8.0)*stardata->common.separation*R_SUN*pow(e1_p2com,-15.0/2.0)*(f1_p - pow(e1_p2com,3.0/2.0)*f2_p*omegaquot2_p);
//		double a1dot_GR_diss_p = (-64.0/5.0)*CONST_GR*m1*m2*(m1+m2)*pow(stardata->common.separation*R_SUN,-3.0)*pow(e1_p2com,-7.0/2.0)*f_GR_adot_p;

		printf("e1dot_third_component = %g; e1dot_tides_star1 = %g; e1dot_tides_star2 = %g; e1dot_GR = %g;\n",e1dot_third_component_p,e1dot_tides_star1_p,e1dot_tides_star2_p,e1dot_GR_diss_p);

		printf("a1dot_tides_star1 = %g; a1dot_tides_star2 = %g; a1dot_GR = %g; ra5 = %g\n",a1dot_tides_star1_p,a1dot_tides_star2_p,a1dot_GR_diss_p,pow(stardata->star[2].radius/stardata->common.separation,8.0));

		printf("--------------------------------------------------------------------------------------------- \n");	

//		PrintFinalStats(cvode_mem);
		printf("Debug mode - press TAB-ENTER to continue to next timestep. ");
		while (1)
		{
			if ('	' == getchar()) break;
		}
		printf("============================================================================================== \n");
	}


	/* Free yev and abstol vectors */
	N_VDestroy_Serial(yev);
	N_VDestroy_Serial(yev_reached);
	N_VDestroy_Serial(abstol);

	/* Free integrator memory */
	CVodeFree(&cvode_mem);

#ifdef IGNORE
	if (output_log == TRUE)
	{
		char cc[70];
		if (change_timestep == FALSE)
		{
			snprintf(cc,40,"TRIPLE %g %g %g",common->eccentricity,common->triple_outer_orbit_eccentricity, \
				common->triple_total_inclination*180.0/PI);
		}
		else
		{
			snprintf(cc,40,"TRIPLE %g %g %g CT %g",common->eccentricity,common->triple_outer_orbit_eccentricity, \
				common->triple_total_inclination*180.0/PI, 100.0*(1.0-timestep_new/timestep_binary_c));
		}

		output_to_logfile(stardata->model.log_fp,
			stardata->model.model_time,
			stardata->star[1].mass,
			stardata->star[2].mass,
			stardata->star[1].stellar_type,
			stardata->star[2].stellar_type,
			stardata->common.separation,
			stardata->common.eccentricity,
			stardata->star[1].radius/stardata->star[1].rol,
			stardata->star[2].radius/stardata->star[2].rol,
			cc,
			stardata);
    	}
#endif

}


static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nSolver final statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

static int fev(realtype t, N_Vector yev, N_Vector ydot, void *data)
{
	
	UserData data_f;
	data_f = (UserData) data;
	struct stardata_t * RESTRICT stardata = data_f->stardata_u;

	Boolean do_print; 
	int iprint;
	
	/*	Constants which appear in the ODE right hand sides	*/
	double a2_f = data_f->a2_initial;				/* Unit: cm */
	double m1_f = data_f->m1;					/* Unit: M_SUN */
	double m2_f = data_f->m2;					/* Unit: M_SUN */
	double m3_f = data_f->m3;					/* Unit: M_SUN */
	double R1_f = data_f->radius_star1;				/* Unit: cm */
	double R2_f = data_f->radius_star2;				/* Unit: cm */
	double amc_star1_f = data_f->apsidal_motion_constant_star1;	/* Apsidal Motion Constant, star 1 */
	double amc_star2_f = data_f->apsidal_motion_constant_star2;	/* Apsidal Motion Constant, star 2 */
	double rg1_f = data_f->gyration_radius_star1;			/* Gyration radius, star 1 */
	double rg2_f = data_f->gyration_radius_star2;			/* Gyration radius, star 2 */
	double kT_star1_f = data_f->kT_star1;				/* Tidal constant, star 1 */
	double kT_star2_f = data_f->kT_star2;				/* Tidal constant, star 2 */
	Boolean RLOF_boolean_f = data_f->RLOF_boolean;

	/*	The ODE variables	*/
//	double e1_f = Ith(yev,1);
//	double e2_f = Ith(yev,2);	
	double x_f = Ith(yev,1);
	double y_f = Ith(yev,2);
	double g1_f = Ith(yev,3);
	double g2_f = Ith(yev,4);
	double a1_f = Ith(yev,5);
//	double Gatot_f = Ith(yev,6);
	double cositot_f = Ith(yev,6);
	double omegaspin1_f = Ith(yev,7);
	double omegaspin2_f = Ith(yev,8);

	double e1_f = 1.0 - pow(10.0,x_f);
	double e2_f = 1.0 - pow(10.0,y_f);

	/*	Some mass combinations	*/
	double q_f = m2_f/m1_f;
	double qinv_f = 1.0/q_f;
	double m_prod_inner_f = m1_f*m2_f;
	double m_tot_inner_f = m1_f+m2_f;
	double m_prod_triple_f = m_prod_inner_f*m3_f;
	double m_tot_triple_f = m_tot_inner_f+m3_f;
	double m_prod_tot_inner_f = m_prod_inner_f*m_tot_inner_f;	/* Used for GR equations */

	/*	Some quantities of interest for tides	*/
	double ra1_f = R1_f/a1_f;
	double ra1_f2 = ra1_f*ra1_f;		/* In the following, the number after the '_f' denotes the power */
	double ra1_f5 = pow(ra1_f,5.0);
	double ra1_f6 = ra1_f*ra1_f5;
	double ra2_f = R2_f/a1_f;
	double ra2_f2 = ra2_f*ra2_f;
	double ra2_f5 = pow(ra2_f,5.0);
	double ra2_f6 = ra2_f*ra2_f5;
	
	double tcqr1_f = q_f*ra1_f6*kT_star1_f;		/* Common combination of quantities in Hut's equations */
	double tcqr2_f = qinv_f*ra2_f6*kT_star2_f;

	double omega_inner_orbit_f = CONST_GM_SQRT*sqrt(m_tot_inner_f*pow(a1_f,-3.0));
	double omegaquot1_f = omegaspin1_f/omega_inner_orbit_f;
	double omegaquot2_f = omegaspin2_f/omega_inner_orbit_f;

	/*	Eccentricity functions	*/
	double e1_f2 = e1_f*e1_f;
	double e1_f4 = e1_f2*e1_f2;
	double e1_f6 = e1_f2*e1_f4;
	double e1_f8 = e1_f4*e1_f4;
	double f1_f = 1.0 + (31.0/2.0)*e1_f2 + (255.0/8.0)*e1_f4 + (185.0/16.0)*e1_f6 + (25.0/64.0)*e1_f8;
	double f2_f = 1.0 + (15.0/2.0)*e1_f2 + (45.0/8.0)*e1_f4 + (5.0/16.0)*e1_f6;
	double f3_f = 1.0 + (15.0/4.0)*e1_f2 + (15.0/8.0)*e1_f4 + (5.0/64.0)*e1_f6;
	double f4_f = 1.0 + (3.0/2.0)*e1_f2 + (1.0/8.0)*e1_f4;
	double f5_f = 1.0 + 3.0*e1_f2 + (3.0/8.0)*e1_f4;
	double f_GR_adot = 1.0 + (73.0/24.0)*e1_f2 + (37.0/96.0)*e1_f4;
	double f_GR_edot = 1.0 + (121.0/304.0)*e1_f2;
	double f_GR_Gtotdot = 1.0 + (7.0/8.0)*e1_f2;
	double e1_f2com = 1.0 - e1_f2;		/* 'com' stands for complement */
	double e2_f2 = e2_f*e2_f;
	double e2_f2com = 1.0 - e2_f2;
	
	/*	Quantities used for triple dynamics	*/
	double L1_f = CONST_ANG_TRIPLE*m_prod_inner_f*sqrt(a1_f/m_tot_inner_f);
	double L2_f = CONST_ANG_TRIPLE*m3_f*m_tot_inner_f*sqrt(a2_f/m_tot_triple_f);
	double Ga1_f = L1_f*sqrt(1.0 - e1_f2);
	double Ga2_f = L2_f*sqrt(1.0 - e2_f2);

	double a1a2quot_f = a1_f/a2_f;
	double C2_f = CONST_GMM*(1.0/16.0)*m_prod_triple_f*pow(m_tot_inner_f,-1.0)*pow(e2_f2com,-1.5)*pow(a1a2quot_f,2.0)*pow(a2_f,-1.0);
	double C3_f = CONST_GMM*(-15.0/16.0)*(1.0/4.0)*m_prod_triple_f*pow(m_tot_inner_f,-2.0)*(m1_f-m2_f)*pow(e2_f2com,-2.5)*pow(a1a2quot_f,3.0)*pow(a2_f,-1.0);
//	printf("C2 C3 %g %g\n",C2_f,C3_f);

//	double cositot_f = (pow(Gatot_f,2.0) - pow(Ga1_f,2.0) - pow(Ga2_f,2.0))/(2.0*Ga1_f*Ga2_f);
	if (cositot_f > 1.0)
	{
		cositot_f = 2.0 - cositot_f;
	}

	if (cositot_f < -1.0)
	{
		cositot_f = -2.0 - cositot_f;
	}

	double cositot_f2 = cositot_f*cositot_f;
	double sinitot_f = sqrt(1.0 - cositot_f2);
	double sinitot_f2 = sinitot_f*sinitot_f;

//	if (isnan(cositot_f)!=0)
//	{
//		printf("cositot = %g \n",cositot_f);
//	}
//	if (isnan(sinitot_f)!=0)
//	{
//		printf("sinitot = %g \n",sinitot_f);
//	}

	double sing1_f = sin(g1_f);
	double sing1_fd = sin(2.0*g1_f);	/* 'd' stands for double angle */
	double sing2_f = sin(g2_f);
	double cosg1_f = cos(g1_f);
	double cosg1_fd = cos(2.0*g1_f);
	double cosg2_f = cos(g2_f);
	
	/*	These quantities appear in the octupole equations	*/
	double B_f = 2.0 + 5.0*e1_f2 - 7.0*e1_f2*cosg1_fd;
	double A_f = 4.0 + 3.0*e1_f2 - (5.0/2.0)*B_f*sinitot_f2;
	double cosphi_f = -cosg1_f*cosg2_f - cositot_f*sing1_f*sing2_f;
	

	/* ================================================================
	* The calculations of the changes of e1, e2, g1, g2, a1, Gatot, 
	* omegaspin1 & omegaspin2
	*/
	
	double factor_triple = 1.0;	/* Set to 0 to disable triple effects */
	double factor_tides_diss = 1.0;	/* Set to 0 to disable dissipative tidal effects in the inner binary */
	double factor_tides_am = 1.0;	/* Set to 0 to disable apsidal motion due to stellar non-sphericity (SDAM) */
	double factor_GR_diss = 1.0;	/* Set to 0 to disable dissipative GR effects, i.e. gravitational wave radiation */
	double factor_GR_am = 1.0;	/* Set to 0 to disable apsidal motion due to GR (GRAM) */

//	if (t*1.0e-6/YEAR_LENGTH_IN_SECONDS > 1245.0)
//	{
//		factor_tides_diss = 0.0;
//		factor_GR_diss = 0.0;
//	}

	/* ----------------------------------------------------------------
	* e1dot: due to third component, tides and GR radiation
	*/ 

	double e1dot_third_component = C2_f*((1.0-e1_f2)/Ga1_f)*(30.0*e1_f*sinitot_f2*sing1_fd) \
		+ C3_f*e2_f*((1.0-e1_f2)/Ga1_f)*(35.0*cosphi_f*sinitot_f2*e1_f2*sing1_fd \
			- 10.0*cositot_f*sinitot_f2*cosg1_f*sing2_f*(1.0-e1_f2) \
			- A_f*(sing1_f*cosg2_f - cositot_f*cosg1_f*sing2_f));
	double e1dot_tides_star1 = -27.0*(1.0+q_f)*tcqr1_f*ra1_f2*e1_f*pow(e1_f2com,-13.0/2.0)*(f3_f \
		- (11.0/18.0)*pow(e1_f2com,3.0/2.0)*f4_f*omegaquot1_f);
	double e1dot_tides_star2 = -27.0*(1.0+qinv_f)*tcqr2_f*ra2_f2*e1_f*pow(e1_f2com,-13.0/2.0)*(f3_f \
		- (11.0/18.0)*pow(e1_f2com,3.0/2.0)*f4_f*omegaquot2_f);
	double e1dot_GR_diss = (-304.0/15.0)*CONST_GR*m_prod_tot_inner_f*e1_f*pow(a1_f,-4.0)*pow(e1_f2com,-5.0/2.0)*f_GR_edot;
//	printf("et %g %g %g %g \n",e1dot_third_component,e1dot_tides_star1,e1dot_tides_star2,e1dot_GR_rad);

//	if ((e1_f <= 1.0e-10) && (100.0*e1dot_third_component < (e1dot_tides_star1 + e1dot_tides_star2 + e1dot_GR_diss)))
//	{
//		printf("ok t %g\n",t*1.0e-6/YEAR_LENGTH_IN_SECONDS);
//	}

	double e1dot_f = factor_triple*e1dot_third_component + factor_tides_diss*(e1dot_tides_star1 + e1dot_tides_star2) \
		+ factor_GR_diss*e1dot_GR_diss;
	Ith(ydot,1) = -1.0*pow(10.0,-1.0*x_f)*e1dot_f/log(10.0);

//	Ith(ydot,1) = factor_triple*e1dot_third_component + factor_tides_diss*(e1dot_tides_star1 + e1dot_tides_star2) \
		+ factor_GR_diss*e1dot_GR_diss;


	/* ----------------------------------------------------------------
	* e2dot: due to third component
	*/ 

	double e2dot_third_component = -C3_f*e1_f*((1.0-e2_f2)/Ga2_f)*(10.0*cositot_f*sinitot_f2*(1.0-e1_f2)*sing1_f*cosg2_f \
		+ A_f*(cosg1_f*sing2_f - cositot_f*sing1_f*cosg2_f));
//	Ith(ydot,2) = factor_triple*e2dot_third_component;

	double e2dot_f = factor_triple*e2dot_third_component;
	Ith(ydot,2) = -1.0*pow(10.0,-1.0*y_f)*e2dot_f/log(10.0);

//	printf("e2d %g \n",e2dot_third_component);
	

	/* ----------------------------------------------------------------
	* g1dot: due to third component, tides and GR
	*/ 

	double g1dot_third_component = 6.0*C2_f*((1.0/Ga1_f)*(4.0*cositot_f2 + (5.0*cosg1_fd - 1.0)*(1.0 - e1_f2 - cositot_f2)) \
			+ (cositot_f/Ga2_f)*(2.0 + e1_f2*(3.0 - 5.0*cosg1_fd))) \
		- C3_f*e2_f*(e1_f*((1.0/Ga2_f) + (cositot_f/Ga1_f))*	(sing1_f*sing2_f*(10.0*(3.0*cositot_f2 \
			- 1.0)*(1.0 - e1_f2) + A_f) - 5.0*B_f*cositot_f*cosphi_f) \
			- ((1.0-e1_f2)/(e1_f*Ga1_f))*(sing1_f*sing2_f*10.0*cositot_f*sinitot_f2*(1.0 - 3.0*e1_f2) \
			+ cosphi_f*(3.0*A_f - 10.0*cositot_f2 + 2.0)));
	double g1dot_tides_star1 = ra1_f5*amc_star1_f*omega_inner_orbit_f*pow(e1_f2com,-2.0)*(15.0*q_f*f4_f*pow(e1_f2com,-3.0) \
		+ (1.0+q_f)*omegaquot1_f*omegaquot1_f);
	double g1dot_tides_star2 = ra2_f5*amc_star2_f*omega_inner_orbit_f*pow(e1_f2com,-2.0)*(15.0*qinv_f*f4_f*pow(e1_f2com,-3.0) \
		+ (1.0+qinv_f)*omegaquot2_f*omegaquot2_f);	
	double g1dot_GR = 3.0*pow(CONST_C_LIGHT_POW2*a1_f*e1_f2com,-1.0)*pow(CONST_GM*m_tot_inner_f/a1_f,3.0/2.0);
//	printf("g1dot %g %g %g %g \n",g1dot_third_component,g1dot_tides_star1,g1dot_tides_star2,g1dot_GR);
	
	Ith(ydot,3) = factor_triple*g1dot_third_component + factor_tides_am*(g1dot_tides_star1 + g1dot_tides_star2) + factor_GR_am*g1dot_GR;


	/* ----------------------------------------------------------------
	* g2dot: due to third component
	*/ 

	double g2dot_third_component = 3.0*C2_f*((2.0*cositot_f/Ga1_f)*(2.0 + e1_f2*(3.0 - 5.0*cosg1_fd)) \
			+ (1.0/Ga2_f)*(4.0 + 6.0*e1_f2 + (5.0*cositot_f2 - 3.0)*(2.0 + e1_f2*(3.0 - 5.0*cosg1_fd)))) \
		+ C3_f*e1_f*(sing1_f*sing2_f*(((4.0*e2_f2 + 1.0)/(e2_f*Ga2_f))*10.0*cositot_f*sinitot_f2*(1.0 - e1_f2) \
			- e2_f*((1.0/Ga1_f) + (cositot_f/Ga2_f))*(A_f + 10.0*(3.0*cositot_f2 - 1.0)*(1.0 - e1_f2))) \
			+ cosphi_f*(5.0*B_f*cositot_f*e2_f*((1.0/Ga1_f) + (cositot_f/Ga2_f)) + ((4.0*e2_f2 + 1.0)/(e2_f*Ga2_f))*A_f));
	Ith(ydot,4) = factor_triple*g2dot_third_component;		
	

	/* ----------------------------------------------------------------
	* a1dot: due to tides and GR radiation
	*/ 

	double a1dot_tides_star1 = -6.0*(1.0+q_f)*tcqr1_f*ra1_f2*a1_f*pow(e1_f2com,-15.0/2.0)*(f1_f \
		- pow(e1_f2com,3.0/2.0)*f2_f*omegaquot1_f);
	double a1dot_tides_star2 = -6.0*(1.0+qinv_f)*tcqr2_f*ra2_f2*a1_f*pow(e1_f2com,-15.0/2.0)*(f1_f \
		- pow(e1_f2com,3.0/2.0)*f2_f*omegaquot2_f);
	double a1dot_GR_diss = (-64.0/5.0)*CONST_GR*m_prod_tot_inner_f*pow(a1_f,-3.0)*pow(e1_f2com,-7.0/2.0)*f_GR_adot;
//	printf("adot %g %g %g\n",a1dot_tides_star1,a1dot_tides_star2,a1dot_GR_rad);

//	if (RLOF_boolean_f == TRUE)
//	{
//		a1dot_tides_star1 = 0;
//		a1dot_tides_star2 = 0;
//	}

	Ith(ydot,5) = factor_tides_diss*(a1dot_tides_star1 + a1dot_tides_star2) + factor_GR_diss*a1dot_GR_diss;


	/* ----------------------------------------------------------------
	* cositotdot: due to triple interaction
	*/ 

	double Ga1dot_third_component = (-1.0)*Ga1_f*e1_f*e1dot_third_component/e1_f2com;
	double Ga2dot_third_component = (-1.0)*Ga2_f*e2_f*e2dot_third_component/e2_f2com;
	double cositotdot_third_component = (-1.0/(Ga1_f*Ga2_f))*(Ga1dot_third_component*(Ga1_f + Ga2_f*cositot_f) + Ga2dot_third_component*(Ga2_f + Ga1_f*cositot_f));
//	if ((cositot_f <= -1.0) && (cositotdot_third_component < 0.0)) cositotdot_third_component *= -1.0;
//	if ((cositot_f >= 1.0) && (cositotdot_third_component > 0.0)) cositotdot_third_component *= -1.0;

	Ith(ydot,6) = factor_triple*cositotdot_third_component;


#ifdef IGNORE
	/* ----------------------------------------------------------------
	* Gatotdot: due to GR radiation and tides in the inner binary
	*/ 

	double Ga1dot_GR_diss = (-32.0/5.0)*CONST_GR*M_SUN*pow(a1_f,-3.0)*m_prod_inner_f*m_prod_inner_f*pow(e1_f2com,-2.0)*pow(CONST_GM*m_tot_inner_f/a1_f,1.0/2.0)*f_GR_Gtotdot;
	double Ga1dot_tides_diss = Ga1_f*(0.5*(a1dot_tides_star1 + a1dot_tides_star2)/a1_f - (e1_f2com/e1_f)*(e1dot_tides_star1 + e1dot_tides_star2)) \
		+ rg1_f*m1_f*M_SUN*R1_f*R1_f*omegaspindot_star1 + rg2_f*m2_f*M_SUN*R2_f*R2_f*omegaspindot_star2;
//	printf("Ga1dot_tides %g\n",Ga1dot_tides_diss);

	Ith(ydot,6) = ((Ga1_f + cositot_f*Ga2_f)/Gatot_f)*(factor_GR_diss*Ga1dot_GR_diss + 0*factor_tides_diss*Ga1dot_tides_diss);
#endif

	/* ----------------------------------------------------------------
	* omegaspindot_star_1: due to tides
	*/ 

	double omegaspindot_star1 = 3.0*tcqr1_f*q_f*pow(rg1_f,-1.0)*omega_inner_orbit_f*pow(e1_f2com,-6.0)*(f2_f \
		- pow(e1_f2com,3.0/2.0)*f5_f*omegaquot1_f);
	Ith(ydot,7) = factor_tides_diss*omegaspindot_star1;


	/* ----------------------------------------------------------------
	* omegaspindot_star_2: due to tides
	*/ 	

	double omegaspindot_star2 = 3.0*tcqr2_f*qinv_f*pow(rg2_f,-1.0)*omega_inner_orbit_f*pow(e1_f2com,-6.0)*(f2_f \
		- pow(e1_f2com,3.0/2.0)*f5_f*omegaquot2_f);
	Ith(ydot,8) = factor_tides_diss*omegaspindot_star2;


	/*	Set do_print to true to output yev and ydot on each call of fev.	*/
	do_print = FALSE;
	if (do_print == TRUE) 
	{
		printf("\n cositot: %g \n",cositot_f);;	
		printf("\n C2, C3: %g %g \n",C2_f, C3_f);	
		for (iprint=1; iprint<5; iprint++) {
			printf("yev %d %g \n",iprint,Ith(yev,iprint));
			printf("ydot %d %g \n",iprint,Ith(ydot,iprint));
		}
//		printf("g1GR %g \n",pararr_f[5]*pow(1.0-e1_f2,-1.0));
		sleep(1.0);

	}



	/* KCTF/KCGWE logging */

//	if (fabs(e1dot_third_component) > 1e-20)
//	{


		/* KCTF */
		

//		if (fabs(e1dot_third_component) < max(fabs(e1dot_tides_star1), fabs(e1dot_tides_star2)))
//		{
		if ((TIDES_CRIT_STAR1 || TIDES_CRIT_STAR2) && (fabs(e1dot_third_component) > 1e-18) && (stardata->common.eccentricity > 0.0) )
		{

			/* KCTF on MS: primary & secondary are MS stars */
			if ( (stardata->star[1].stellar_type <= 1) && (stardata->star[2].stellar_type <= 1) && (stardata->common.triple_KCTF[5] == 0) ) 
			{
//				printf("TC t %g\n",stardata->model.model_time);
				stardata->common.triple_KCTF_data_before[5][0] = stardata->model.model_time;
				stardata->common.triple_KCTF_data_before[5][1] = stardata->common.separation/AURSUN; /*AU*/
				stardata->common.triple_KCTF_data_before[5][2] = stardata->common.triple_outer_semi_major_axis; /*AU*/
				stardata->common.triple_KCTF_data_before[5][3] = log10(1.0-stardata->common.eccentricity);
				stardata->common.triple_KCTF_data_before[5][4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);				
				stardata->common.triple_KCTF_data_before[5][5] = stardata->common.triple_total_inclination*180.0/M_PI; /*deg*/
				stardata->common.triple_KCTF_data_before[5][6] = stardata->star[1].mass;
				stardata->common.triple_KCTF_data_before[5][7] = stardata->star[2].mass;
				stardata->common.triple_KCTF_data_before[5][8] = stardata->common.triple_outer_mass;
				stardata->common.triple_KCTF_data_before_st[5][1] = stardata->star[1].stellar_type;
				stardata->common.triple_KCTF_data_before_st[5][2] = stardata->star[2].stellar_type;
				stardata->common.triple_KCTF[5] ++;
			}
	
//			if ( (log10(fabs(stardata->common.triple_a1dot_tides_star1)) < -2) && (log10(fabs(stardata->common.triple_a1dot_tides_star2)) < -2) && (stardata->common.triple_KCTF[5] == 1) )
			if (stardata->common.triple_KCTF[5] == 1)

			{
				printf("id %d KCTF MS t %g \n",stardata->common.triple_thread_id,stardata->model.model_time);
				stardata->common.triple_KCTF_data_after[5][0] = stardata->model.model_time;
				stardata->common.triple_KCTF_data_after[5][1] = stardata->common.separation/AURSUN; /*AU*/
				stardata->common.triple_KCTF_data_after[5][2] = stardata->common.triple_outer_semi_major_axis; /*AU*/
				stardata->common.triple_KCTF_data_after[5][3] = log10(1.0-stardata->common.eccentricity);
				stardata->common.triple_KCTF_data_after[5][4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);				
				stardata->common.triple_KCTF_data_after[5][5] = stardata->common.triple_total_inclination*180.0/M_PI; /*deg*/
				stardata->common.triple_KCTF_data_after[5][6] = stardata->star[1].mass;
				stardata->common.triple_KCTF_data_after[5][7] = stardata->star[2].mass;
				stardata->common.triple_KCTF_data_after[5][8] = stardata->common.triple_outer_mass;
				stardata->common.triple_KCTF_data_after_st[5][1] = stardata->star[1].stellar_type;
				stardata->common.triple_KCTF_data_after_st[5][2] = stardata->star[2].stellar_type;
				stardata->common.triple_KCTF[5] ++;

				char buff_KCTF[50] = "";
				sprintf(buff_KCTF,"KCTF_data%d.txt",stardata->common.triple_thread_id);
				FILE *fp_KCTF;
				fp_KCTF = fopen(buff_KCTF,"a");
				fprintf(fp_KCTF,"5\t%g\t%g\t%g\t%g\t%g\t\%g\t%g\t%g\t%g\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\n",stardata->common.triple_KCTF_data_before[5][0],stardata->common.triple_KCTF_data_before[5][1],stardata->common.triple_KCTF_data_before[5][2],stardata->common.triple_KCTF_data_before[5][3],stardata->common.triple_KCTF_data_before[5][4],stardata->common.triple_KCTF_data_before[5][5],stardata->common.triple_KCTF_data_before[5][6],stardata->common.triple_KCTF_data_before[5][7],stardata->common.triple_KCTF_data_before[5][8],stardata->common.triple_KCTF_data_before_st[5][1],stardata->common.triple_KCTF_data_before_st[5][2],stardata->common.triple_KCTF_data_after[5][0],stardata->common.triple_KCTF_data_after[5][1],stardata->common.triple_KCTF_data_after[5][2],stardata->common.triple_KCTF_data_after[5][3],stardata->common.triple_KCTF_data_after[5][4],stardata->common.triple_KCTF_data_after[5][5],stardata->common.triple_KCTF_data_after[5][6],stardata->common.triple_KCTF_data_after[5][7],stardata->common.triple_KCTF_data_after[5][8],stardata->common.triple_KCTF_data_after_st[5][1],stardata->common.triple_KCTF_data_after_st[5][2]);
				fclose(fp_KCTF);
			}

			/* KCTF: compact object + MS */
			if ( (stardata->star[1].stellar_type >= 10) && (stardata->star[1].stellar_type <= 13) && (stardata->star[2].stellar_type <= 1) && (stardata->common.triple_KCTF[6] == 0) ) 
			{
//				printf("TC t %g\n",stardata->model.model_time);
				stardata->common.triple_KCTF_data_before[6][0] = stardata->model.model_time;
				stardata->common.triple_KCTF_data_before[6][1] = stardata->common.separation/AURSUN; /*AU*/
				stardata->common.triple_KCTF_data_before[6][2] = stardata->common.triple_outer_semi_major_axis; /*AU*/
				stardata->common.triple_KCTF_data_before[6][3] = log10(1.0-stardata->common.eccentricity);
				stardata->common.triple_KCTF_data_before[6][4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);				
				stardata->common.triple_KCTF_data_before[6][5] = stardata->common.triple_total_inclination*180.0/M_PI; /*deg*/
				stardata->common.triple_KCTF_data_before[6][6] = stardata->star[1].mass;
				stardata->common.triple_KCTF_data_before[6][7] = stardata->star[2].mass;
				stardata->common.triple_KCTF_data_before[6][8] = stardata->common.triple_outer_mass;
				stardata->common.triple_KCTF_data_before_st[6][1] = stardata->star[1].stellar_type;
				stardata->common.triple_KCTF_data_before_st[6][2] = stardata->star[2].stellar_type;
				stardata->common.triple_KCTF[6] ++;
			}
//			if ( (log10(fabs(stardata->common.triple_a1dot_tides_star2)) < -5) && (stardata->common.triple_KCTF[5] == 1) )
			if (stardata->common.triple_KCTF[6] == 1)

			{
				printf("id %d KCTF Compact object + MS t %g \n",stardata->common.triple_thread_id,stardata->model.model_time);
				stardata->common.triple_KCTF_data_after[6][0] = stardata->model.model_time;
				stardata->common.triple_KCTF_data_after[6][1] = stardata->common.separation/AURSUN; /*AU*/
				stardata->common.triple_KCTF_data_after[6][2] = stardata->common.triple_outer_semi_major_axis; /*AU*/
				stardata->common.triple_KCTF_data_after[6][3] = log10(1.0-stardata->common.eccentricity);
				stardata->common.triple_KCTF_data_after[6][4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);				
				stardata->common.triple_KCTF_data_after[6][5] = stardata->common.triple_total_inclination*180.0/M_PI; /*deg*/
				stardata->common.triple_KCTF_data_after[6][6] = stardata->star[1].mass;
				stardata->common.triple_KCTF_data_after[6][7] = stardata->star[2].mass;
				stardata->common.triple_KCTF_data_after[6][8] = stardata->common.triple_outer_mass;
				stardata->common.triple_KCTF_data_after_st[6][1] = stardata->star[1].stellar_type;
				stardata->common.triple_KCTF_data_after_st[6][2] = stardata->star[2].stellar_type;
				stardata->common.triple_KCTF[6] ++;
	
				char buff_KCTF[50] = "";
				sprintf(buff_KCTF,"KCTF_data%d.txt",stardata->common.triple_thread_id);
				FILE *fp_KCTF;
				fp_KCTF = fopen(buff_KCTF,"a");
				fprintf(fp_KCTF,"6\t%g\t%g\t%g\t%g\t%g\t\%g\t%g\t%g\t%g\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\n",stardata->common.triple_KCTF_data_before[6][0],stardata->common.triple_KCTF_data_before[6][1],stardata->common.triple_KCTF_data_before[6][2],stardata->common.triple_KCTF_data_before[6][3],stardata->common.triple_KCTF_data_before[6][4],stardata->common.triple_KCTF_data_before[6][5],stardata->common.triple_KCTF_data_before[6][6],stardata->common.triple_KCTF_data_before[6][7],stardata->common.triple_KCTF_data_before[6][8],stardata->common.triple_KCTF_data_before_st[6][1],stardata->common.triple_KCTF_data_before_st[6][2],stardata->common.triple_KCTF_data_after[6][0],stardata->common.triple_KCTF_data_after[6][1],stardata->common.triple_KCTF_data_after[6][2],stardata->common.triple_KCTF_data_after[6][3],stardata->common.triple_KCTF_data_after[6][4],stardata->common.triple_KCTF_data_after[6][5],stardata->common.triple_KCTF_data_after[6][6],stardata->common.triple_KCTF_data_after[6][7],stardata->common.triple_KCTF_data_after[6][8],stardata->common.triple_KCTF_data_after_st[6][1],stardata->common.triple_KCTF_data_after_st[6][2]);
				fclose(fp_KCTF);
			}
		}

		/* KCGWE */
//		if ((fabs(e1dot_third_component) < fabs(e1dot_GR_diss)) && (x_f < -4))
		if ((GWE_CRIT) && (fabs(e1dot_third_component) > 1e-18) )
		{
			if (stardata->common.triple_KCTF[7] == 0) 
			{
				printf("KCGWE t %g\n",stardata->model.model_time);
				stardata->common.triple_KCTF_data_before[7][0] = stardata->model.model_time;
				stardata->common.triple_KCTF_data_before[7][1] = stardata->common.separation/AURSUN; /*AU*/
				stardata->common.triple_KCTF_data_before[7][2] = stardata->common.triple_outer_semi_major_axis; /*AU*/
				stardata->common.triple_KCTF_data_before[7][3] = log10(1.0-stardata->common.eccentricity);
				stardata->common.triple_KCTF_data_before[7][4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);				
				stardata->common.triple_KCTF_data_before[7][5] = stardata->common.triple_total_inclination*180.0/M_PI; /*deg*/
				stardata->common.triple_KCTF_data_before[7][6] = stardata->star[1].mass;
				stardata->common.triple_KCTF_data_before[7][7] = stardata->star[2].mass;
				stardata->common.triple_KCTF_data_before[7][8] = stardata->common.triple_outer_mass;
				stardata->common.triple_KCTF_data_before_st[7][1] = stardata->star[1].stellar_type;
				stardata->common.triple_KCTF_data_before_st[7][2] = stardata->star[2].stellar_type;
				stardata->common.triple_KCTF[7] ++;
			}
		
			if ((fabs(a1dot_GR_diss) < 1e-5) && (stardata->common.triple_KCTF[7] == 1))
			{
				printf("id %d KCGWE t %g \n",stardata->common.triple_thread_id,stardata->model.model_time);
				stardata->common.triple_KCTF_data_after[7][0] = stardata->model.model_time;
				stardata->common.triple_KCTF_data_after[7][1] = stardata->common.separation/AURSUN; /*AU*/
				stardata->common.triple_KCTF_data_after[7][2] = stardata->common.triple_outer_semi_major_axis; /*AU*/
				stardata->common.triple_KCTF_data_after[7][3] = log10(1.0-stardata->common.eccentricity);
				stardata->common.triple_KCTF_data_after[7][4] = log10(1.0-stardata->common.triple_outer_orbit_eccentricity);				
				stardata->common.triple_KCTF_data_after[7][5] = stardata->common.triple_total_inclination*180.0/M_PI; /*deg*/
				stardata->common.triple_KCTF_data_after[7][6] = stardata->star[1].mass;
				stardata->common.triple_KCTF_data_after[7][7] = stardata->star[2].mass;
				stardata->common.triple_KCTF_data_after[7][8] = stardata->common.triple_outer_mass;
				stardata->common.triple_KCTF_data_after_st[7][1] = stardata->star[1].stellar_type;
				stardata->common.triple_KCTF_data_after_st[7][2] = stardata->star[2].stellar_type;
				stardata->common.triple_KCTF[7] ++;

				char buff_KCTF[50] = "";
				sprintf(buff_KCTF,"KCTF_data%d.txt",stardata->common.triple_thread_id);
				FILE *fp_KCTF;
				fp_KCTF = fopen(buff_KCTF,"a");
				fprintf(fp_KCTF,"7\t%g\t%g\t%g\t%g\t%g\t\%g\t%g\t%g\t%g\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\n",stardata->common.triple_KCTF_data_before[7][0],stardata->common.triple_KCTF_data_before[7][1],stardata->common.triple_KCTF_data_before[7][2],stardata->common.triple_KCTF_data_before[7][3],stardata->common.triple_KCTF_data_before[7][4],stardata->common.triple_KCTF_data_before[7][5],stardata->common.triple_KCTF_data_before[7][6],stardata->common.triple_KCTF_data_before[7][7],stardata->common.triple_KCTF_data_before[7][8],stardata->common.triple_KCTF_data_before_st[7][1],stardata->common.triple_KCTF_data_before_st[7][2],stardata->common.triple_KCTF_data_after[7][0],stardata->common.triple_KCTF_data_after[7][1],stardata->common.triple_KCTF_data_after[7][2],stardata->common.triple_KCTF_data_after[7][3],stardata->common.triple_KCTF_data_after[7][4],stardata->common.triple_KCTF_data_after[7][5],stardata->common.triple_KCTF_data_after[7][6],stardata->common.triple_KCTF_data_after[7][7],stardata->common.triple_KCTF_data_after[7][8],stardata->common.triple_KCTF_data_after_st[7][1],stardata->common.triple_KCTF_data_after_st[7][2]);
				fclose(fp_KCTF);
			}
		}
	



	return(0);
}


static int froot(realtype t, N_Vector yev, realtype *gout, void *data)
{
	UserData data_f;
	data_f = (UserData) data;

//	gout[0] = Ith(yev,5) - (1.0 - (A1_CHANGE_TOL/100.0))*data_f->a1_initial;
//	gout[1] = Ith(yev,5) - (1.0 + (A1_CHANGE_TOL/100.0))*data_f->a1_initial;
//	gout[2] = Ith(yev,7) - (1.0 - (OMSPIN1_CHANGE_TOL/100.0))*data_f->omega_spin_star1_initial;
//	gout[3] = Ith(yev,7) - (1.0 + (OMSPIN1_CHANGE_TOL/100.0))*data_f->omega_spin_star1_initial;
//	gout[4] = Ith(yev,8) - (1.0 - (OMSPIN2_CHANGE_TOL/100.0))*data_f->omega_spin_star2_initial;
//	gout[5] = Ith(yev,8) - (1.0 + (OMSPIN2_CHANGE_TOL/100.0))*data_f->omega_spin_star2_initial;

	double e1_t = 1.0 - pow(10.0,Ith(yev,1)); /* Current e1 */
	double e2_t = 1.0 - pow(10.0,Ith(yev,2)); /* Current e2 */

	gout[0] = Ith(yev,5)*(1.0 - e1_t)/(data_f->radius_star1 + data_f->radius_star2) - 1.0;

	double G1changefraction;
	G1changefraction = 0.02;
	if (data_f->ST1 <= 1)
	{
		G1changefraction = 0.2;
	}
//	printf("ST1 = %d; G1cf = %g\n",data_f->ST1,G1changefraction);

//	gout[1] = fabs( sqrt( (Ith(yev,5)*(1.0 - Ith(yev,1)*Ith(yev,1))) / (data_f->a1_initial*(1.0 - data_f->e1_initial*data_f->e1_initial)) ) - 1.0) - 0.02;
//	gout[1] = fabs( sqrt( (Ith(yev,5)*(1.0 - data_f->e1_initial*data_f->e1_initial)) / (data_f->a1_initial*(1.0 - data_f->e1_initial*data_f->e1_initial)) ) - 1.0) - 0.02;
	gout[1] = fabs( sqrt( (Ith(yev,5)*(1.0 - e1_t*e1_t)) / (data_f->a1_initial*(1.0 - data_f->e1_initial*data_f->e1_initial)) ) - 1.0) - G1changefraction;

	return(0);
}

static void ehfun(int error_code, const char *module, const char *function, char *msg, void *eh_data)
{
	UserData data_f;
	data_f = (UserData) eh_data;
	struct stardata_t * RESTRICT stardata = data_f->stardata_u;

	int thread_id = data_f->thread_id;
	
	/*	id errcode time a1 x ST1 ST2 m1 m2 m3 RL1 RL2 e1dot_STD e1dot_TF1 e1dot_TF2 e1dot_GR a1dot_TF1 a1dot_TF2 a1dot_GR	*/

	FILE *errorfile;
	errorfile = fopen("CVODE_errors.txt","a");
	fprintf(errorfile,"%d %d %g %g %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g \n",thread_id,error_code,data_f->t_start/(1e6*YEAR_LENGTH_IN_SECONDS),data_f->a1_initial/(AURSUN*R_SUN),data_f->a2_initial/(AURSUN*R_SUN),log10(1.0 - data_f->e1_initial),data_f->ST1,data_f->ST2,stardata->star[1].mass,stardata->star[2].mass,stardata->common.triple_outer_mass,stardata->star[1].radius*pow(stardata->star[1].rol,-1.0),stardata->star[2].radius*pow(stardata->star[2].rol,-1.0),stardata->common.triple_e1dot_third_component,stardata->common.triple_e1dot_tides_star1,stardata->common.triple_e1dot_tides_star2,stardata->common.triple_e1dot_GR_diss,stardata->common.triple_a1dot_tides_star1,stardata->common.triple_a1dot_tides_star2,stardata->common.triple_a1dot_GR_diss,msg);
	fclose(errorfile);

	printf("error %s %s %s\n",module,function,msg);
	printf("error code; thread_id %d %d\n",error_code,thread_id);
}
