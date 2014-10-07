/*	Worker code for SecularTriple, a secular triple gravitational dynamics code taking into account Newtonian, 1PN and 2.5PN terms	*/
/*	The relevant ODEs are solved consistently for each user supplied timestep using CVODE (Cohen & Hindmarsh 1996)	*/

#include "main_code.h"

//#define NEQ		    10				        /* number of ODE equations */
#define RTOL		RCONST(1.0e-10)			/* scalar relative tolerance - acceptable value is 1.0e-10 as determined by trial and error */
#define ATOL1		RCONST(1.0e-8)
#define ATOL2		RCONST(1.0e-8)
#define ATOL3		RCONST(1.0e-8)
#define ATOL4		RCONST(1.0e-8)
#define ATOL5		RCONST(1.0e1)
#define ATOL6		RCONST(1.0e1)
#define ATOL7		RCONST(1.0e-8)
#define ATOL8		RCONST(1.0e-8)	
#define ATOL9		RCONST(1.0e-8)				
#define ATOL10		RCONST(1.0e-8)				
#define TIMESTEP0	RCONST(1.0)				/* initial internal ODE timestep */
#define MAXNUMSTEPS 	5e8					/* maximum number of internal steps */
#define MAXTIME		RCONST(13.7e9*365.25*24.0*3600.0)	/* maximum integration time */
#define MAXNUMCONVFAIL	20					/* maximum number of convergence failures */

/*	interface parameters	*/
int equations_of_motion_specification = 0;
bool check_for_dynamical_stability;
bool check_for_inner_collision,check_for_outer_collision;
bool check_for_inner_RLOF,check_for_outer_RLOF;
bool include_quadrupole_terms,include_octupole_terms;
bool include_1PN_inner_terms,include_1PN_outer_terms,include_1PN_inner_outer_terms,include_25PN_inner_terms,include_25PN_outer_terms;
bool include_inner_tidal_terms,include_outer_tidal_terms;
bool include_inner_wind_terms,include_outer_wind_terms;
bool include_wind_spin_coupling_terms;
bool include_inner_RLOF_terms,include_outer_RLOF_terms;
bool include_linear_mass_change,include_linear_radius_change;
double relative_tolerance = 1.0e-10;

int get_equations_of_motion_specification(int *equations_of_motion_specification_t)
{
    *equations_of_motion_specification_t = equations_of_motion_specification;
    return 0;
}
int set_equations_of_motion_specification(int equations_of_motion_specification_t)
{
    equations_of_motion_specification = equations_of_motion_specification_t;
    return 0;
}
int get_relative_tolerance(double *relative_tolerance_t)
{
    *relative_tolerance_t = relative_tolerance;
    return 0;
}
int set_relative_tolerance(double relative_tolerance_t)
{
    relative_tolerance = relative_tolerance_t;
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
int get_include_inner_wind_terms(int *value){
    *value = include_inner_wind_terms ? 1 : 0;
    return 0;
}
int set_include_inner_wind_terms(int value){
    include_inner_wind_terms = value == 1;
    return 0;
}
int get_include_outer_wind_terms(int *value){
    *value = include_outer_wind_terms ? 1 : 0;
    return 0;
}
int set_include_outer_wind_terms(int value){
    include_outer_wind_terms = value == 1;
    return 0;
}
int get_include_wind_spin_coupling_terms(int *value){
    *value = include_wind_spin_coupling_terms ? 1 : 0;
    return 0;
}
int set_include_wind_spin_coupling_terms(int value){
    include_wind_spin_coupling_terms = value == 1;
    return 0;
}
int get_include_inner_RLOF_terms(int *value){
    *value = include_inner_RLOF_terms ? 1 : 0;
    return 0;
}
int set_include_inner_RLOF_terms(int value){
    include_inner_RLOF_terms = value == 1;
    return 0;
}
int get_include_outer_RLOF_terms(int *value){
    *value = include_outer_RLOF_terms ? 1 : 0;
    return 0;
}
int set_include_outer_RLOF_terms(int value){
    include_outer_RLOF_terms = value == 1;
    return 0;
}
int get_include_linear_mass_change(int *value){
    *value = include_linear_mass_change ? 1 : 0;
    return 0;
}
int set_include_linear_mass_change(int value){
    include_linear_mass_change = value == 1;
    return 0;
}
int get_include_linear_radius_change(int *value){
    *value = include_linear_radius_change ? 1 : 0;
    return 0;
}
int set_include_linear_radius_change(int value){
    include_linear_radius_change = value == 1;
    return 0;
}

int evolve(
    int stellar_type1, int stellar_type2, int stellar_type3,
    double m1, double m2, double m3,
    double m1_convective_envelope, double m2_convective_envelope, double m3_convective_envelope,
    double R1, double R2, double R3,
    double R1_convective_envelope, double R2_convective_envelope, double R3_convective_envelope,    
    double luminosity_star1, double luminosity_star2, double luminosity_star3,
    double spin_angular_frequency1, double spin_angular_frequency2, double spin_angular_frequency3,
    double AMC_star1, double AMC_star2, double AMC_star3,
    double gyration_radius_star1, double gyration_radius_star2, double gyration_radius_star3,
//    double k_div_T_tides_star1, double k_div_T_tides_star2, double k_div_T_tides_star3,
    double a_in, double a_out,
    double e_in, double e_out,
    double INCL_in, double INCL_out, double AP_in, double AP_out, double LAN_in, double LAN_out,
    bool star1_is_donor, bool star2_is_donor, bool star3_is_donor,
    double wind_mass_loss_rate_star1, double wind_mass_loss_rate_star2, double wind_mass_loss_rate_star3,
    double time_derivative_of_radius_star1, double time_derivative_of_radius_star2, double time_derivative_of_radius_star3,
    double inner_mass_transfer_rate, double outer_mass_transfer_rate,
    double inner_accretion_efficiency_wind_child1_to_child2, double inner_accretion_efficiency_wind_child2_to_child1,
    double outer_accretion_efficiency_wind_child1_to_child2,double outer_accretion_efficiency_wind_child2_to_child1,
    double inner_accretion_efficiency_mass_transfer, double outer_accretion_efficiency_mass_transfer,
    double inner_specific_AM_loss_mass_transfer, double outer_specific_AM_loss_mass_transfer,    
    double t, double dt,
    double * m1_output, double * m2_output, double * m3_output,
    double * R1_output, double * R2_output, double * R3_output,
    double * spin_angular_frequency1_output, double * spin_angular_frequency2_output, double * spin_angular_frequency3_output,
    double * a_in_output, double * a_out_output,
    double * e_in_output, double * e_out_output,
    double *INCL_in_output, double *INCL_out_output, double *INCL_in_out_output, double * AP_in_output, double * AP_out_output, double *LAN_in_output, double *LAN_out_output,
    double * t_output,
    int * output_flag, int * error_flag)
{
    double tiny_double = 1.0e-14;
    if (e_in==0.0) { e_in = tiny_double; }
    if (e_out==0.0) {e_out = tiny_double; }
    if (INCL_in==0.0) { INCL_in = tiny_double; }
    if (INCL_out==0.0) { INCL_out = tiny_double; }
    if (AP_in==0.0) { AP_in = tiny_double; }
    if (AP_out==0.0) { AP_out = tiny_double; }
    if (LAN_in==0.0) { LAN_in = tiny_double; }
    if (LAN_out==0.0) { LAN_out = tiny_double; }

    if (e_in<=tiny_double) { e_in = tiny_double; }

    /*********************************************************************
     * quantities that are assumed to be constant during the integration *
     *********************************************************************/
	UserData data;
	data = NULL;
	data = (UserData) malloc(sizeof *data);
    
    data->dt = dt; // the global time-step
    
    data->stellar_type1 = stellar_type1;
    data->stellar_type2 = stellar_type2;
    data->stellar_type3 = stellar_type3;
	data->m1 = m1;
	data->m2 = m2;
	data->m3 = m3;
    data->m1_convective_envelope = m1_convective_envelope;
    data->m2_convective_envelope = m2_convective_envelope;
    data->m3_convective_envelope = m3_convective_envelope;
	data->R1 = R1;
	data->R2 = R2;
    data->R3 = R3;
    data->R1_convective_envelope = R1_convective_envelope;
    data->R2_convective_envelope = R2_convective_envelope;
    data->R3_convective_envelope = R3_convective_envelope;

    data->include_quadrupole_terms = include_quadrupole_terms;
    data->include_octupole_terms = include_octupole_terms;    
    data->include_1PN_inner_terms = include_1PN_inner_terms;
    data->include_1PN_outer_terms = include_1PN_outer_terms;
    data->include_1PN_inner_outer_terms = include_1PN_inner_outer_terms;
    data->include_25PN_inner_terms = include_25PN_inner_terms;
    data->include_25PN_outer_terms = include_25PN_outer_terms;

    data->include_inner_tidal_terms = include_inner_tidal_terms;
    data->include_outer_tidal_terms = include_outer_tidal_terms;
    data->AMC_star1 = AMC_star1; // Apsidal Motion Constant
    data->AMC_star2 = AMC_star2; // Apsidal Motion Constant
    data->AMC_star3 = AMC_star3; // Apsidal Motion Constant
    data->gyration_radius_star1 = gyration_radius_star1; // gyration radius (NOT squared)
    data->gyration_radius_star2 = gyration_radius_star2; // gyration radius (NOT squared)
    data->gyration_radius_star3 = gyration_radius_star3; // gyration radius (NOT squared)
    data->luminosity_star1 = luminosity_star1;
    data->luminosity_star2 = luminosity_star2;
    data->luminosity_star3 = luminosity_star3;
//    data->k_div_T_tides_star1 = k_div_T_tides_star1; // tidal dissipation constant
//    data->k_div_T_tides_star2 = k_div_T_tides_star2; // tidal dissipation constant
//    data->k_div_T_tides_star3 = k_div_T_tides_star3; // tidal dissipation constant

    data->include_linear_mass_change = include_linear_mass_change;
    data->include_linear_radius_change = include_linear_radius_change;

    data->include_inner_wind_terms = include_inner_wind_terms;
    data->include_outer_wind_terms = include_outer_wind_terms;
    data->include_inner_RLOF_terms = include_inner_RLOF_terms;
    data->include_outer_RLOF_terms = include_outer_RLOF_terms;
    data->include_wind_spin_coupling_terms = include_wind_spin_coupling_terms;
    data->star1_is_donor = star1_is_donor;
    data->star2_is_donor = star2_is_donor;
    data->star3_is_donor = star3_is_donor;
    data->wind_mass_loss_rate_star1 = wind_mass_loss_rate_star1;
    data->wind_mass_loss_rate_star2 = wind_mass_loss_rate_star2;
    data->wind_mass_loss_rate_star3 = wind_mass_loss_rate_star3;
    data->time_derivative_of_radius_star1 = time_derivative_of_radius_star1;
    data->time_derivative_of_radius_star2 = time_derivative_of_radius_star2;
    data->time_derivative_of_radius_star3 = time_derivative_of_radius_star3;
    data->inner_mass_transfer_rate = inner_mass_transfer_rate;
    data->outer_mass_transfer_rate = outer_mass_transfer_rate;
    data->inner_accretion_efficiency_wind_child1_to_child2 = inner_accretion_efficiency_wind_child1_to_child2;
    data->inner_accretion_efficiency_wind_child2_to_child1 = inner_accretion_efficiency_wind_child2_to_child1;
    data->outer_accretion_efficiency_wind_child1_to_child2 = outer_accretion_efficiency_wind_child1_to_child2;
    data->outer_accretion_efficiency_wind_child2_to_child1 = outer_accretion_efficiency_wind_child2_to_child1;
    data->inner_accretion_efficiency_mass_transfer = inner_accretion_efficiency_mass_transfer;
    data->outer_accretion_efficiency_mass_transfer = outer_accretion_efficiency_mass_transfer;
    data->inner_specific_AM_loss_mass_transfer = inner_specific_AM_loss_mass_transfer;
    data->outer_specific_AM_loss_mass_transfer = outer_specific_AM_loss_mass_transfer;

    data->check_for_dynamical_stability = check_for_dynamical_stability;
    data->check_for_inner_collision = check_for_inner_collision;
    data->check_for_outer_collision = check_for_outer_collision;
    data->check_for_inner_RLOF = check_for_inner_RLOF;
    data->check_for_outer_RLOF = check_for_outer_RLOF;

//	double L1 = m1*m2*sqrt(CONST_G*a1/(m1+m2));
//	double L2 = (m1+m2)*m3*sqrt(CONST_G*a2/(m1+m2+m3));
//	double e1_2com = 1.0 - e1*e1;
//	double e2_2com = 1.0 - e2*e2;
//	double Ga1 = L1*sqrt(e1_2com);
//	double Ga2 = L2*sqrt(e2_2com);
//	double Gatot =  sqrt(pow(Ga1,2.0) + pow(Ga2,2.0) + 2.0*Ga1*Ga2*cos(i12));	

	N_Vector yev, yev_out, abstol;
	void *cvode_mem;
	int flag,flag_s;

	yev = yev_out = abstol = NULL;
	cvode_mem = NULL;

    int NEQ;
    double e_in_vec[3], e_out_vec[3], h_in_vec[3], h_out_vec[3], q_in_vec[3], q_out_vec[3];    
    double cos_INCL_in = cos(INCL_in);
    double cos_INCL_out = cos(INCL_out);
    double sin_INCL_in = sin(INCL_in);
    double sin_INCL_out = sin(INCL_out);
    double cos_INCL_in_out = sin_INCL_in*sin_INCL_out*cos(LAN_in-LAN_out) + cos_INCL_in*cos_INCL_out;
    double INCL_in_out = acos(cos_INCL_in_out);

    if (equations_of_motion_specification==0)
    {
        /********************************************************
         * solve equations of motion based on delaunay elements *
         * ******************************************************/

        NEQ = 10;
        yev = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev, "N_VNew_Serial", 0)) return;
        yev_out = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev_reached, "N_VNew_Serial", 0)) return;
        abstol = N_VNew_Serial(NEQ); 
//	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return;

        Ith(yev,1) = log10(1.0 - e_in);
        Ith(yev,2) = log10(1.0 - e_out);
        Ith(yev,3) = AP_in;		/* Inner orbit argument of periastron - unit: rad */
        Ith(yev,4) = AP_out;    /* Outer orbit argument of periastron - unit: rad */
        Ith(yev,5) = a_in;		/* Inner orbit semi-major axis - unit: cm */
        Ith(yev,6) = a_out;		/* Outer orbit semi-major axis - unit: cm */
        Ith(yev,7) = cos_INCL_in_out;
        Ith(yev,8) = spin_angular_frequency1;	/* Spin angular frequency of star 1 - unit: rad s^-1 */
        Ith(yev,9) = spin_angular_frequency2;	/* Spin angular frequency of star 2 - unit: rad s^-1 */
        Ith(yev,10) = spin_angular_frequency3;	/* Spin angular frequency of star 3 - unit: rad s^-1 */

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

        cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
//        flag = CVodeInit(cvode_mem, fev_delaunay, t, yev);       
        double start_time = 0.0;
        flag = CVodeInit(cvode_mem, fev_delaunay, start_time, yev);
        
    }
    else
    {
        /********************************************************************
         * solve equations of motion based on triad of unit vectors (q,h,e) *
         * ******************************************************************/        
         
        NEQ = 18;
        yev = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev, "N_VNew_Serial", 0)) return;
        yev_out = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev_reached, "N_VNew_Serial", 0)) return;
        abstol = N_VNew_Serial(NEQ); 
//	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return;         

        compute_triad_vectors_from_orbital_elements(m1,m2,m3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out,
            e_in_vec,e_out_vec,h_in_vec,h_out_vec,q_in_vec,q_out_vec);
        double x_in_vec[3], x_out_vec[3];
        int i=0;
        for (i=0; i<3; i++)
        {
            x_in_vec[i] = log10(1.0-e_in)*e_in_vec[i]/e_in;
            x_out_vec[i] = log10(1.0-e_out)*e_out_vec[i]/e_out;
            
            Ith(yev,i+1) = x_in_vec[i];
            Ith(yev,i+4) = x_out_vec[i];
            Ith(yev,i+7) = h_in_vec[i];
            Ith(yev,i+10) = h_out_vec[i];
            Ith(yev,i+13) = 0.0; // reserved for stellar spin vectors
            Ith(yev,i+16) = 0.0;
        }
        
        FILE *fp;
        fp = fopen("/data2/amuse-svn/output_stream.txt","w");
        fprintf(fp,"test %g %g %g ",e_in_vec[2],h_in_vec[2],q_in_vec[2]);
        fclose(fp);

        cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        flag = CVodeInit(cvode_mem, fev_triad, t, yev);            

    }


    /**************************/
    /***    CVode setup     ***/
    /**************************/

//	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
//	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return;
	
//	flag = CVodeInit(cvode_mem, fev, t, yev);
//	if (check_flag(&flag, "CVodeInit", 1)) return;

//	flag = CVodeSetErrHandlerFn(cvode_mem, ehfun, eh_data);
//	if (check_flag(&flag, "CVodeSetErrHandlerFn", 1)) return;

//	flag = CVodeSVtolerances(cvode_mem, RTOL, abstol);
    flag = CVodeSVtolerances(cvode_mem, relative_tolerance, abstol);
//	if (check_flag(&flag, "CVodeSVtolerances", 1)) return;
	
	flag = CVDense(cvode_mem, NEQ);
//	if (check_flag(&flag, "CVDense", 1)) return;
  
	flag = CVodeSetInitStep(cvode_mem, TIMESTEP0);
//	if (check_flag(&flag, "CVodeSetInitStep", 1)) return;
		
	flag = CVodeSetMaxNumSteps(cvode_mem, MAXNUMSTEPS);
//	if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return;

//	flag = CVodeSetMinStep(cvode_mem, 0.1);
//	if (check_flag(&flag, "CVodeSetMinStep", 1)) return;

	flag = CVodeSetMaxHnilWarns(cvode_mem, 1);
//	if (check_flag(&flag, "CVodeSetMaxHnilWarns", 1)) return;
			
//	flag = CVodeSetStopTime(cvode_mem, MAXTIME);
//	if (check_flag(&flag, "CVodeSetStopTime", 1)) return(1);

	flag = CVodeSetMaxConvFails(cvode_mem, MAXNUMCONVFAIL);
//	if (check_flag(&flag, "CVodeSetMaxConvFails", 1)) return;
			
	flag = CVodeSetUserData(cvode_mem, data);
//	if (check_flag(&flag, "CVodeSetUsetData", 1)) return;

	int nroot = 6;
	int rootsfound[nroot];
	flag = CVodeRootInit(cvode_mem, nroot, froot_delaunay);
//	if (check_flag(&flag, "CVodeRootInit", 1)) return;	

	double t_end = t + dt;
	double t_end_cvode;
//	int output_flag = 0;
//	int error_flag = 0;

//	flag_s = CVode(cvode_mem, t_end, yev_out, &t_end_cvode, CV_NORMAL);	
	flag_s = CVode(cvode_mem, dt, yev_out, &t_end_cvode, CV_NORMAL);	
    
	if (flag_s == CV_SUCCESS)
	{
		*output_flag = 0;
		*error_flag = 0;
	}
	else if (flag_s == CV_ROOT_RETURN)
	{
		CVodeGetRootInfo(cvode_mem,rootsfound);
		if (rootsfound[0] == 1 || rootsfound[0] == -1) // dynamical instability
		{
			*output_flag = 1;
		}
		if (rootsfound[1] == 1 || rootsfound[1] == -1) // collision in inner binary
		{
			*output_flag = 2;
		}
		if (rootsfound[2] == 1 || rootsfound[2] == -1) // collision in outer binary -- see ODE_system.c for what this means
		{
			*output_flag = 3;
		}    
		if (rootsfound[3] == 1 || rootsfound[3] == -1) // RLOF star1 -- think about increasing case
		{
			*output_flag = 4;
		}    
		if (rootsfound[4] == 1 || rootsfound[4] == -1) // RLOF star2 -- think about increasing case
		{
			*output_flag = 5;
		}    
		if (rootsfound[5] == 1 || rootsfound[5] == -1) // RLOF star3 -- think about increasing case
		{
			*output_flag = 6;
		}    
		*error_flag = 0;
	}
	else
		*output_flag = 99;
		*error_flag = flag_s;


    /**********
     * output *
     * ********/
     
	*m1_output = m1;
	*m2_output = m2;
	*m3_output = m3;
	*R1_output = R1;
	*R2_output = R2;
    *R3_output = R3;
    
	*t_output = t_end_cvode;
//	*output_flag = output_flag;
//	*error_flag = error_flag;

    double x_output,y_output;
//      double L1_out,L2_out,G1_out,G2_out,G12_out;
    double cos_INCL_in_output;
    double cos_INCL_out_output;
    double sin_INCL_in_output;
    double sin_INCL_out_output;
    double cos_INCL_in_out_output;
//    double INCL12 = acos(cos_INCL12);

    if (equations_of_motion_specification==0)
    {
        x_output = Ith(yev_out,1);
        y_output = Ith(yev_out,2);
        
        *e_in_output = 1.0 - pow(10,x_output);
        *e_out_output = 1.0 - pow(10,y_output);
        *AP_in_output = Ith(yev_out,3);
        *AP_out_output = Ith(yev_out,4);
        *a_in_output = Ith(yev_out,5);
        *a_out_output = Ith(yev_out,6);        

        cos_INCL_in_out_output = Ith(yev_out,7);

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
        
        *spin_angular_frequency1_output = Ith(yev_out,8);
        *spin_angular_frequency2_output = Ith(yev_out,9);
        *spin_angular_frequency3_output = Ith(yev_out,10);
        
//        L1_out = m1*m2*sqrt(CONST_G*(*a1_out)/(m1+m2));
//        L2_out = (m1+m2)*m3*sqrt(CONST_G*(*a2_out)/(m1+m2+m3));
//        G1_out = L1_out*sqrt(1.0 - (*e1_out)*(*e1_out));
//        G2_out = L2_out*sqrt(1.0 - (*e2_out)*(*e2_out));
//        G12_out = sqrt(G1_out*G1_out + G2_out*G2_out + 2.0*G1_out*G2_out*theta_out);
//        cos_i1_out = (G12_out*G12_out + G1_out*G1_out - G2_out*G2_out)/(2.0*G12_out*G1_out);
//        cos_i2_out = (G12_out*G12_out - G1_out*G1_out + G2_out*G2_out)/(2.0*G12_out*G2_out);
        
//        *i1_out = acos(cos_i1_out);
//        *i2_out = acos(cos_i2_out);
    }
    else
    {
#ifdef IGNORE        
        compute_orbital_elements_from_triad_vectors(m1,m2,m3,
            e1_vec,e2_vec,h1_vec,h2_vec,q1_vec,q2_vec,
            a1_out,a2_out,e1_out,e2_out,
            INCL1_out,INCL2_out,AP1_out,AP2_out,LAN1_out,LAN2_out);
            
        cos_INCL1_out = cos(*INCL1_out);
        cos_INCL2_out = cos(*INCL2_out);
        sin_INCL1_out = sin(*INCL1_out);
        sin_INCL2_out = sin(*INCL2_out);
        cos_INCL12_out = sin_INCL1_out*sin_INCL2_out*cos(*LAN1_out-*LAN2_out) + cos_INCL1_out*cos_INCL2_out;
        *INCL12_out = acos(cos_INCL12_out);

        FILE *fp;
        fp = fopen("/data2/amuse-svn/output_stream3.txt","w");
        fprintf(fp,"test %g %g %g %g ",*a2_out/a2,*e2_out/e2,cos(*INCL1_out),cos(INCL1));
        fclose(fp);
        
        if (1==0)
        {
        *a1_out = a1;
        *a2_out = a2;
        *e1_out = e1;
        *e2_out = e2;
        *INCL1_out = INCL1;
        *INCL2_out = INCL2;
        *INCL12_out = INCL12;
        *AP1_out = AP1;
        *AP2_out = AP2;
        *LAN1_out = LAN1;
        *LAN2_out = LAN2;
        }
#endif
    }

	return 0;
}

