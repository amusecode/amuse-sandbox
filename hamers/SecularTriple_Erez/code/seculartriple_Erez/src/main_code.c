/*	Worker code for SecularTriple, a secular triple gravitational dynamics code taking into account Newtonian, 1PN and 2.5PN terms	*/
/*	The relevant ODEs are solved consistently for each user supplied timestep using CVODE (Cohen & Hindmarsh 1996)	*/

#include "main_code.h"

//#define NEQ		    10				        /* number of ODE equations */
#define RTOL		RCONST(1.0e-10)			/* scalar relative tolerance - acceptable value is 1.0e-10 as determined by trial and error */
#define ATOL1		RCONST(1.0e-8)
#define ATOL2		RCONST(1.0e-8)
#define ATOL3		RCONST(1.0e-8)
#define ATOL4		RCONST(1.0e-8)
#define ATOL5		RCONST(1.0e-8)
#define ATOL6		RCONST(1.0e-8)
#define ATOL7		RCONST(1.0e-8)
#define ATOL8		RCONST(1.0e-8)	
#define ATOL9		RCONST(1.0e-8)				
#define ATOL10		RCONST(1.0e-2)				
#define ATOL11		RCONST(1.0e-2)				
#define ATOL12		RCONST(1.0e-2)				
#define INITIAL_ODE_TIME_STEP	RCONST(1.0e-10)				/* initial internal ODE timestep */
#define MAXNUMSTEPS 	5e8					/* maximum number of internal steps */
#define MAXTIME		RCONST(13.7e9*365.25*24.0*3600.0)	/* maximum integration time */
#define MAXNUMCONVFAIL	20					/* maximum number of convergence failures */
#define MAXERRTESTFAILS 20

/*	interface parameters	*/
int equations_of_motion_specification = 0;
bool check_for_dynamical_stability;
bool check_for_inner_collision,check_for_outer_collision;
bool check_for_inner_RLOF,check_for_outer_RLOF;
bool include_quadrupole_terms,include_octupole_terms;
bool include_1PN_inner_terms,include_1PN_outer_terms,include_1PN_inner_outer_terms,include_25PN_inner_terms,include_25PN_outer_terms;
bool include_inner_tidal_terms,include_outer_tidal_terms;

bool ignore_tertiary;
bool verbose = TRUE;
double relative_tolerance = 1.0e-10;
double input_precision = 1.0e-5;
double threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = 1.0e-12;
double threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero = 1.0e-8;
int linear_solver = 0;

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
    double t, double global_time_step,
    double * m1_output, double * m2_output, double * m3_output,
    double * R1_output, double * R2_output, double * R3_output,
    double * tidal_E1_dot_output, double * tidal_E2_dot_output, double * tidal_E3_dot_output,
    double * spin_angular_frequency1_output, double * spin_angular_frequency2_output, double * spin_angular_frequency3_output,
    double * a_in_output, double * a_out_output,
    double * e_in_output, double * e_out_output,
    double *INCL_in_output, double *INCL_out_output, double *INCL_in_out_output, double * AP_in_output, double * AP_out_output, double *LAN_in_output, double *LAN_out_output,
    double * t_output,
    int * CVODE_flag, int * root_finding_flag
)
{
    double tiny_double = input_precision;
    if (e_in<=tiny_double) { e_in = tiny_double; }
    if (e_out<=tiny_double) { e_out = tiny_double; }    
    if (INCL_in<=tiny_double) { INCL_in = tiny_double; }
    if (INCL_out<=tiny_double) { INCL_out = tiny_double; }

    if (INCL_in > M_PI-tiny_double)
    {
        INCL_in = M_PI-tiny_double;
    }
    if (INCL_out > M_PI-tiny_double)
    {
        INCL_out = M_PI-tiny_double;
    }

//    if (spin_angular_frequency1<=tiny_double) { spin_angular_frequency1 = tiny_double; }
//    if (spin_angular_frequency2<=tiny_double) { spin_angular_frequency2 = tiny_double; }
//    if (spin_angular_frequency3<=tiny_double) { spin_angular_frequency3 = tiny_double; }
    
    /*********************************************************************
     * ODE parameters *
     *********************************************************************/
	UserData data;
	data = NULL;
	data = (UserData) malloc(sizeof *data);
    
    data->global_time_step = global_time_step;
    
	data->m1 = m1;
	data->m2 = m2;
	data->m3 = m3;
	data->R1 = R1;
	data->R2 = R2;
    data->R3 = R3;

    data->include_quadrupole_terms = include_quadrupole_terms;
    data->include_octupole_terms = include_octupole_terms;    
    data->include_1PN_inner_terms = include_1PN_inner_terms;
    data->include_1PN_outer_terms = include_1PN_outer_terms;
    data->include_1PN_inner_outer_terms = include_1PN_inner_outer_terms;
    data->include_25PN_inner_terms = include_25PN_inner_terms;
    data->include_25PN_outer_terms = include_25PN_outer_terms;
    data->ignore_tertiary = ignore_tertiary;

    data->threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero;
    data->include_inner_tidal_terms = include_inner_tidal_terms;
    data->include_outer_tidal_terms = include_outer_tidal_terms;
    data->AMC_star1 = AMC_star1; // Apsidal Motion Constant
    data->AMC_star2 = AMC_star2; // Apsidal Motion Constant
    data->AMC_star3 = AMC_star3; // Apsidal Motion Constant
    data->gyration_radius_star1 = gyration_radius_star1; // gyration radius (NOT squared)
    data->gyration_radius_star2 = gyration_radius_star2; // gyration radius (NOT squared)
    data->gyration_radius_star3 = gyration_radius_star3; // gyration radius (NOT squared)

    data->k_div_T_tides_star1 = k_div_T_tides_star1; // tidal dissipation constant
    data->k_div_T_tides_star2 = k_div_T_tides_star2; // tidal dissipation constant
    data->k_div_T_tides_star3 = k_div_T_tides_star3; // tidal dissipation constant

    data->check_for_dynamical_stability = check_for_dynamical_stability;
    data->check_for_inner_collision = check_for_inner_collision;
    data->check_for_outer_collision = check_for_outer_collision;
    data->check_for_inner_RLOF = check_for_inner_RLOF;
    data->check_for_outer_RLOF = check_for_outer_RLOF;

	N_Vector yev, yev_out, abstol;
	void *cvode_mem;
	int flag,flag_s;

	yev = yev_out = abstol = NULL;
	cvode_mem = NULL;

    int NEQ;

    double INCL_in_out = INCL_in;
    double cos_INCL_in_out = cos(INCL_in_out);

    double start_time = 0.0;
    
    if (equations_of_motion_specification==0)
    {
        /********************************************************
         * solve equations of motion based on delaunay elements *
         * ******************************************************/

        NEQ = 12;
        yev = N_VNew_Serial(NEQ);
        if (check_flag((void *)yev, "N_VNew_Serial", 0)) return 1;
        yev_out = N_VNew_Serial(NEQ);
        if (check_flag((void *)yev_out, "N_VNew_Serial", 0)) return 1;
        abstol = N_VNew_Serial(NEQ); 
    	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return 1;

        Ith(yev,1) = log10(1.0 - e_in);
        Ith(yev,2) = log10(1.0 - e_out);
        Ith(yev,3) = AP_in;		/* Inner orbit argument of periastron - unit: rad */
        Ith(yev,4) = AP_out;    /* Outer orbit argument of periastron - unit: rad */
        Ith(yev,5) = LAN_in;		/* Inner orbit longitude of the ascending node - unit: rad */
        Ith(yev,6) = LAN_out;		/* Outer orbit longitude of the ascending node - unit: rad */
        Ith(yev,7) = a_in;		/* Inner orbit semi-major axis - unit: AU */
        Ith(yev,8) = a_out;		/* Outer orbit semi-major axis - unit: AU */
        Ith(yev,9) = cos_INCL_in_out; /* cosine of the relative inclination */
        Ith(yev,10) = spin_angular_frequency1;	/* Spin angular frequency of star 1 - unit: rad Myr^-1 */
        Ith(yev,11) = spin_angular_frequency2;	/* Spin angular frequency of star 2 - unit: rad Myr^-1 */
        Ith(yev,12) = spin_angular_frequency3;	/* Spin angular frequency of star 3 - unit: rad Myr^-1 */

        Ith(abstol,1) = ATOL1;   
        Ith(abstol,2) = ATOL2;
        Ith(abstol,3) = ATOL3;
        Ith(abstol,4) = ATOL4;
        Ith(abstol,5) = ATOL5;   
        Ith(abstol,6) = ATOL6;
        Ith(abstol,7) = ATOL7;
        Ith(abstol,8) = ATOL8;
        Ith(abstol,9) = ATOL9;
        Ith(abstol,10) = ATOL10;
        Ith(abstol,11) = ATOL11;
        Ith(abstol,12) = ATOL12;

        /* use Backward Differentiation Formulas (BDF)
            scheme in conjunction with Newton iteration --
            these choices are recommended for stiff ODEs
            in the CVODE manual                          
        */
        cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        flag = CVodeInit(cvode_mem, fev_delaunay, start_time, yev);
        
    }
    else
    {
        exit(-1);
    }


    /**************************/
    /***    CVode setup     ***/
    /**************************/

	flag = CVodeSetErrHandlerFn(cvode_mem, error_handling_function, data);
	if (check_flag(&flag, "CVodeSetErrHandlerFn", 1)) return 1;

    flag = CVodeSVtolerances(cvode_mem, relative_tolerance, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return 1;

    /****** specify linear solver ******/
    if (linear_solver == 0) // default solver
    {
        flag = CVDense(cvode_mem, NEQ);
    	if (check_flag(&flag, "CVDense", 1)) return 1;
    }
    if (linear_solver == 1)
    {
        flag = CVDiag(cvode_mem);
        if (check_flag(&flag, "CVDiag", 1)) return 1;
    }
    if (linear_solver == 2) // placeholder for CVBand; not yet implemented
    {
        //flag = CVBand(cvode_mem);
        //if (check_flag(&flag, "CVBand", 1)) return 1;
    }
    if (linear_solver == 3)
    {
        flag = CVSpgmr(cvode_mem, PREC_NONE, CVSPILS_MAXL);
        if (check_flag(&flag, "CVSpgmr", 1)) return 1;

    }
    if (linear_solver == 4)
    {
        flag = CVSpbcg(cvode_mem, PREC_NONE, CVSPILS_MAXL);
        if (check_flag(&flag, "CVSpbcg", 1)) return 1;
    }
    if (linear_solver == 5)
    {
        flag = CVSptfqmr(cvode_mem, PREC_NONE, CVSPILS_MAXL);
        if (check_flag(&flag, "CVSptfqmr", 1)) return 1;
    }

	flag = CVodeSetInitStep(cvode_mem, INITIAL_ODE_TIME_STEP);
	if (check_flag(&flag, "CVodeSetInitStep", 1)) return 1;
		
	flag = CVodeSetMaxNumSteps(cvode_mem, MAXNUMSTEPS);
	if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return 1;

	flag = CVodeSetMinStep(cvode_mem, 1.0e-15);
	if (check_flag(&flag, "CVodeSetMinStep", 1)) return 1;

	flag = CVodeSetMaxHnilWarns(cvode_mem, 1);
	if (check_flag(&flag, "CVodeSetMaxHnilWarns", 1)) return 1;
			
//	flag = CVodeSetStopTime(cvode_mem, MAXTIME);
//	if (check_flag(&flag, "CVodeSetStopTime", 1)) return 1;

	flag = CVodeSetMaxConvFails(cvode_mem, MAXNUMCONVFAIL);
	if (check_flag(&flag, "CVodeSetMaxConvFails", 1)) return 1;

	flag = CVodeSetMaxErrTestFails(cvode_mem, MAXERRTESTFAILS);
	if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return 1;

	flag = CVodeSetUserData(cvode_mem, data);
	if (check_flag(&flag, "CVodeSetUserData", 1)) return 1;

	flag = CVodeSetStabLimDet(cvode_mem, TRUE);
	if (check_flag(&flag, "CVodeSetStabLimDet", 1)) return 1;

	flag = CVodeSetMaxOrd(cvode_mem, 5);
	if (check_flag(&flag, "CVodeSetMaxOrd", 1)) return 1;

	int nroot = 6;
	int rootsfound[nroot];
	flag = CVodeRootInit(cvode_mem, nroot, froot_delaunay);
	if (check_flag(&flag, "CVodeRootInit", 1)) return 1;	

	double t_end = start_time + global_time_step;
	double t_end_cvode;

	flag_s = CVode(cvode_mem, t_end, yev_out, &t_end_cvode, CV_NORMAL);	
    
	if (flag_s == CV_SUCCESS) // successfully integrated for global time-step
	{
		*CVODE_flag = 0;
        *root_finding_flag = 0;        
	}
	else if (flag_s == CV_ROOT_RETURN) // root was found
	{
		CVodeGetRootInfo(cvode_mem,rootsfound);
		if (rootsfound[0] == 1 || rootsfound[0] == -1) // dynamical instability
		{
			*root_finding_flag = 1;
		}
		if (rootsfound[1] == 1 || rootsfound[1] == -1) // collision in inner binary
		{
			*root_finding_flag = 2;
		}
		if (rootsfound[2] == 1 || rootsfound[2] == -1) // collision in outer binary -- see ODE_system.c for what this means
		{
			*root_finding_flag = 3;
		}    
        if (rootsfound[3] == 1)
        {
            *root_finding_flag = 4; // star 1 now fills its Roche Lobe: R-R_L was negative, has become positive
        }
        if (rootsfound[3] == -1)
        {
            *root_finding_flag = -4; // star 1 no longer fills its Roche Lobe: R-R_L was positive, has become negative
        }
        if (rootsfound[4] == 1)
        {
            *root_finding_flag = 5; // star 2 now fills its Roche Lobe: R-R_L was negative, has become positive
        }
        if (rootsfound[4] == -1)
        {
            *root_finding_flag = -5; // star 2 no longer fills its Roche Lobe: R-R_L was positive, has become negative
        }
        if (rootsfound[5] == 1)
        {
            *root_finding_flag = 6; // star 3 now fills its Roche Lobe: R-R_L was negative, has become positive
        }
        if (rootsfound[5] == -1)
        {
            *root_finding_flag = -6; // star 3 no longer fills its Roche Lobe: R-R_L was positive, has become negative
        }
		*CVODE_flag = CV_ROOT_RETURN; // indicates root was found
	}
	else if (flag_s == CV_WARNING) // integration successfull, but warnings occured
    {
		*CVODE_flag = CV_WARNING;
		*root_finding_flag = flag_s;
    }
    else if (flag_s < 0) // fatal error(s) occurred
    {
		*CVODE_flag = flag_s;
		*root_finding_flag = 0.0;
    }
    else
    {
        printf("unknown CVODE output code %d\n",flag_s);
        exit(-1);
    }

    /**********
     * output *
     * ********/
     
	*m1_output = m1;
	*m2_output = m2;
	*m3_output = m3;
	*R1_output = R1;
	*R2_output = R2;
    *R3_output = R3;

    *tidal_E1_dot_output = data->tidal_E1_dot;
    *tidal_E2_dot_output = data->tidal_E2_dot;
    *tidal_E3_dot_output = data->tidal_E3_dot;
    
	*t_output = t_end_cvode;

    double x_output,y_output;
    double cos_INCL_in_output;
    double cos_INCL_out_output;
    double sin_INCL_in_output;
    double sin_INCL_out_output;
    double cos_INCL_in_out_output;

    if (equations_of_motion_specification==0)
    {
        x_output = Ith(yev_out,1);
        y_output = Ith(yev_out,2);
        
        *e_in_output = 1.0 - pow(10,x_output);
        *e_out_output = 1.0 - pow(10,y_output);
        *AP_in_output = Ith(yev_out,3);
        *AP_out_output = Ith(yev_out,4);
        *LAN_in_output = Ith(yev_out,5);
        *LAN_out_output = Ith(yev_out,6);
        *a_in_output = Ith(yev_out,7);
        *a_out_output = Ith(yev_out,8);        

        if (*e_in_output<=tiny_double) { *e_in_output = tiny_double; }
        if (*e_out_output<=tiny_double) { *e_out_output = tiny_double; }    

        cos_INCL_in_out_output = Ith(yev_out,9);

        if (cos_INCL_in_out_output > 1.0)
        {
            cos_INCL_in_out_output = 2.0 - cos_INCL_in_out_output;
        }
    
        if (cos_INCL_in_out_output < -1.0)
        {
            cos_INCL_in_out_output = -2.0 - cos_INCL_in_out_output;
        }

        *INCL_in_output = acos(cos_INCL_in_out_output);
        *INCL_out_output = 0.0;
        *INCL_in_out_output = acos(cos_INCL_in_out_output);
        
        *spin_angular_frequency1_output = Ith(yev_out,10);
        *spin_angular_frequency2_output = Ith(yev_out,11);
        *spin_angular_frequency3_output = Ith(yev_out,12);
    }
    else
    {
        exit(-1);
    }

    N_VDestroy_Serial(yev);
    N_VDestroy_Serial(yev_out);
    N_VDestroy_Serial(abstol);
    CVodeFree(&cvode_mem);

	return 0;
}

/* function to check ODE solver-related function return values */
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

void error_handling_function(int error_code, const char *module, const char *function, char *message, void *data_f)
{
    printf("error_handling_function error %d\n",error_code);
    printf("module %s\n",module);
    printf("message %s\n",message);
    UserData data;
	data = (UserData) data_f;
//    data->stop_after_error_bool = TRUE;
}

int get_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero(double *value)
{
    *value = threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero;
    return 0;
}
int set_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero(double value)
{
    threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = value;
    return 0;
}
int get_equations_of_motion_specification(int *value)
{
    *value = equations_of_motion_specification;
    return 0;
}
int set_equations_of_motion_specification(int value)
{
    equations_of_motion_specification = value;
    return 0;
}
int get_relative_tolerance(double *value)
{
    *value = relative_tolerance;
    return 0;
}
int set_relative_tolerance(double value)
{
    relative_tolerance = value;
    return 0;
}
int get_input_precision(double *value)
{
    *value = relative_tolerance;
    return 0;
}
int set_input_precision(double value)
{
    input_precision = value;
    return 0;
}
int get_linear_solver(int *value)
{
    *value = linear_solver;
    return 0;
}
int set_linear_solver(int value)
{
    linear_solver = value;
    return 0;
}
int get_check_for_dynamical_stability(int *value){
    *value = check_for_dynamical_stability ? 1 : 0;
    return 0;
}
int set_check_for_dynamical_stability(int value){
    check_for_dynamical_stability = value == 1;
    return 0;
}
int get_check_for_inner_collision(int *value){
    *value = check_for_inner_collision ? 1 : 0;
    return 0;
}
int set_check_for_inner_collision(int value){
    check_for_inner_collision = value == 1;
    return 0;
}
int get_check_for_outer_collision(int *value){
    *value = check_for_outer_collision ? 1 : 0;
    return 0;
}
int set_check_for_outer_collision(int value){
    check_for_outer_collision = value == 1;
    return 0;
}
int get_check_for_inner_RLOF(int *value){
    *value = check_for_inner_RLOF ? 1 : 0;
    return 0;
}
int set_check_for_inner_RLOF(int value){
    check_for_inner_RLOF = value == 1;
    return 0;
}
int get_check_for_outer_RLOF(int *value){
    *value = check_for_outer_RLOF ? 1 : 0;
    return 0;
}
int set_check_for_outer_RLOF(int value){
    check_for_outer_RLOF = value == 1;
    return 0;
}
int get_include_quadrupole_terms(int *value){
    *value = include_quadrupole_terms ? 1 : 0;
    return 0;
}
int set_include_quadrupole_terms(int value){
    include_quadrupole_terms = value == 1;
    return 0;
}
int get_include_octupole_terms(int *value){
    *value = include_octupole_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_terms(int value){
    include_octupole_terms = value == 1;
    return 0;
}
int get_include_1PN_inner_terms(int *value){
    *value = include_1PN_inner_terms ? 1 : 0;
    return 0;
}
int set_include_1PN_inner_terms(int value){
    include_1PN_inner_terms = value == 1;
    return 0;
}
int get_include_1PN_outer_terms(int *value){
    *value = include_1PN_outer_terms ? 1 : 0;
    return 0;
}
int set_include_1PN_outer_terms(int value){
    include_1PN_outer_terms = value == 1;
    return 0;
}
int get_include_1PN_inner_outer_terms(int *value){
    *value = include_1PN_inner_outer_terms ? 1 : 0;
    return 0;
}
int set_include_1PN_inner_outer_terms(int value){
    include_1PN_inner_outer_terms = value == 1;
    return 0;
}
int get_include_25PN_inner_terms(int *value){
    *value = include_25PN_inner_terms ? 1 : 0;
    return 0;
}
int set_include_25PN_inner_terms(int value){
    include_25PN_inner_terms = value == 1;
    return 0;
}
int get_include_25PN_outer_terms(int *value){
    *value = include_25PN_outer_terms ? 1 : 0;
    return 0;
}
int set_include_25PN_outer_terms(int value){
    include_25PN_outer_terms = value == 1;
    return 0;
}
int get_include_inner_tidal_terms(int *value){
    *value = include_inner_tidal_terms ? 1 : 0;
    return 0;
}
int set_include_inner_tidal_terms(int value){
    include_inner_tidal_terms = value == 1;
    return 0;
}
int get_include_outer_tidal_terms(int *value){
    *value = include_outer_tidal_terms ? 1 : 0;
    return 0;
}
int set_include_outer_tidal_terms(int value){
    include_outer_tidal_terms = value == 1;
    return 0;
}
int get_ignore_tertiary(int *value){
    *value = ignore_tertiary ? 1 : 0;
    return 0;
}
int set_ignore_tertiary(int value){
    ignore_tertiary = value == 1;
    return 0;
}
int get_verbose(int *value){
    *value = verbose ? 1 : 0;
    return 0;
}
int set_verbose(int value){
    verbose = value == 1;
    return 0;
}