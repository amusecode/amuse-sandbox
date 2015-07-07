/*	Worker code for SecularTriple, a secular triple gravitational dynamics code taking into account Newtonian, 1PN and 2.5PN terms	*/
/*  Also included: effects of tidal friction, wind mass loss & mass transfer */
/*	The relevant ODEs are solved consistently for each user supplied timestep using CVODE (Cohen & Hindmarsh 1996)	*/

#include "main_code.h"

int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data_f)
{
    
	UserData data;
	data = (UserData) data_f;
	
    /**********************
     * extract parameters *
     **********************/
    double global_time_step = data->global_time_step; // the global time-step
    
	double m1 = data->m1; // masses -- at the END of the global time-step
	double m2 = data->m2;				
	double m3 = data->m3;		
    double R1 = data->R1; // radii -- at the END of the global time-step
    double R2 = data->R2;
    double R3 = data->R3;

    double AMC_star1 = data->AMC_star1; // Apsidal Motion Constant
    double AMC_star2 = data->AMC_star2; // Apsidal Motion Constant
    double AMC_star3 = data->AMC_star3; // Apsidal Motion Constant  
    double gyration_radius_star1 = data->gyration_radius_star1; // gyration radius (NOT squared)     
    double gyration_radius_star2 = data->gyration_radius_star2; // gyration radius (NOT squared)     
    double gyration_radius_star3 = data->gyration_radius_star3; // gyration radius (NOT squared)             

    double k_div_T_tides_star1 = data->k_div_T_tides_star1; // AMC divided by tidal dissipation time-scale
    double k_div_T_tides_star2 = data->k_div_T_tides_star2;    
    double k_div_T_tides_star3 = data->k_div_T_tides_star3;

    bool include_quadrupole_terms = data->include_quadrupole_terms;
    bool include_octupole_terms = data->include_octupole_terms;
    bool include_1PN_inner_terms = data->include_1PN_inner_terms;
    bool include_1PN_outer_terms = data->include_1PN_outer_terms;
    bool include_1PN_inner_outer_terms = data->include_1PN_inner_outer_terms;
    bool include_25PN_inner_terms = data->include_25PN_inner_terms;
    bool include_25PN_outer_terms = data->include_25PN_outer_terms;
    bool include_inner_tidal_terms = data->include_inner_tidal_terms;
    bool include_outer_tidal_terms = data->include_outer_tidal_terms;
    bool ignore_tertiary = data->ignore_tertiary;
    
    double threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = data->threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero;

	/*	the ODE variables	
     *  g: argument of pericentre
     *  h: longitude of the ascending nodes
     */
     
	double x = Ith(yev,1); // log_10(1-e_in)
	double y = Ith(yev,2); // log_10(1-e_out)
	double g_in = Ith(yev,3);
	double g_out = Ith(yev,4);
    double h_in = Ith(yev,5);
    double h_out = Ith(yev,6);
	double a_in = Ith(yev,7);
	double a_out = Ith(yev,8);
	double cositot = Ith(yev,9);
	double spin_angular_frequency1 = Ith(yev,10);
	double spin_angular_frequency2 = Ith(yev,11);
	double spin_angular_frequency3 = Ith(yev,12);

	double e_in = 1.0 - pow(10.0,x);
	double e_out = 1.0 - pow(10.0,y);

    /* in the case that the tertiary is not to be taken into account (ignore_tertiary == TRUE),
     * several quantities should still be set to some arbitrary, nonzero value in order to avoid nans at various instances
     * this does not affect the values of the quantities outside of this function,
     * nor the dots of these quantities computed below */
    if (ignore_tertiary == TRUE)
    {
        m3 = 1.0;
        a_out = 1.0e10;
        e_out = 0.1;
    }
    
    /********************************
     * preamble: derived quantities *
     * *****************************/
     
	/*	mass quantities	*/
	double m1_div_m2 = m1/m2;
	double m2_div_m1 = 1.0/m1_div_m2;    
    double m1_plus_m2 = m1+m2;
    double m1_plus_m2_plus_m3 = m1_plus_m2+m3;
    double m1_times_m2 = m1*m2;
    double m1_times_m2_times_m3 = m1_times_m2*m3;
    double m1_plus_m2_div_m3 = m1_plus_m2/m3;
    double m1_times_m1 = m1*m1;
    double m2_times_m2 = m2*m2;
        
	/*	eccentricity quantities	*/
	double e_in_p2 = e_in*e_in;
	double l_in_p2 = 1.0 - e_in_p2;
    double l_in = sqrt(l_in_p2);
    double l_in_p3 = l_in_p2*l_in;
    double l_in_p4 = l_in_p3*l_in;
    double l_in_p5 = l_in_p4*l_in;
    double l_in_p6 = l_in_p5*l_in;
    double l_in_p7 = l_in_p6*l_in;
    
	double e_out_p2 = e_out*e_out;
	double l_out_p2 = 1.0 - e_out_p2;
    double l_out = sqrt(l_out_p2);
    double l_out_p3 = l_out_p2*l_out;
    double l_out_p4 = l_out_p3*l_out;
    double l_out_p5 = l_out_p4*l_out;    
    double l_out_p6 = l_out_p5*l_out;
    double l_out_p7 = l_out_p6*l_out;

	/*	triple secular gravitational dynamics quantities */
    /* 2000ApJ...535..385F */
	double L_in = m1_times_m2*sqrt(CONST_G*a_in/(m1_plus_m2));
	double L_out = (m1_plus_m2)*m3*sqrt(CONST_G*a_out/(m1_plus_m2_plus_m3));
	double G_in = L_in*sqrt(l_in_p2);
	double G_out = L_out*sqrt(l_out_p2);
    double G_tot = sqrt( G_in*G_in + G_out*G_out + 2.0*G_in*G_out*cositot );

	double a_in_div_a_out = a_in/a_out;
    double C2,C3;
    if (include_quadrupole_terms == FALSE)
    {
        C2 = 0.0;
    }
    else
    {
        C2 = CONST_G*c_1div16*(m1_times_m2_times_m3/m1_plus_m2)*pow(l_out,-3.0)*a_in_div_a_out*a_in_div_a_out/a_out;
    }
    if (include_octupole_terms == FALSE)
    {
        C3 = 0.0;
    }
    else
    {
        C3 = -CONST_G*c_15div16*c_1div4*(m1_times_m2_times_m3/(m1_plus_m2*m1_plus_m2))*(m1-m2)*pow(l_out,-5.0)*a_in_div_a_out*a_in_div_a_out*a_in_div_a_out/a_out;
    }

	if (cositot > 1.0)
	{
		cositot = 2.0 - cositot;
	}

	if (cositot < -1.0)
	{
		cositot = -2.0 - cositot;
	}

	double cositot_p2 = cositot*cositot;
	double sinitot = sqrt(1.0 - cositot_p2); // NOTE: 0 < itot < PI, so sinitot > 0 always
	double sinitot_p2 = sinitot*sinitot;

	double sin_g_in = sin(g_in);
	double sin_2g_in = sin(2.0*g_in);
	double sin_g_out = sin(g_out);
	double cos_g_in = cos(g_in);
	double cos_2g_in = cos(2.0*g_in);
	double cos_g_out = cos(g_out);
	
	/*	required for octupole-order terms	*/
	double B = 2.0 + 5.0*e_in_p2 - 7.0*e_in_p2*cos_2g_in;
	double A = 4.0 + 3.0*e_in_p2 - c_5div2*B*sinitot_p2;
	double cosphi = -cos_g_in*cos_g_out - cositot*sin_g_in*sin_g_out;

    /* PN quantities */
    double f_25PN_e_in = 0.0, f_25PN_a_in = 0.0;
    if (include_25PN_inner_terms == TRUE)
    {
        f_25PN_e_in = f_25PN_e(e_in_p2);
        f_25PN_a_in = f_25PN_a(e_in_p2);
    }
    double f_25PN_e_out = 0.0, f_25PN_a_out = 0.0;
    if (include_25PN_outer_terms == TRUE)
    {
        f_25PN_e_out = f_25PN_e(e_out_p2);
        f_25PN_a_out = f_25PN_a(e_out_p2);
    }
    double f_1PN_in_out_m1_m2 = 0.0, f_1PN_in_out_L_in_tilde = 0.0, f_1PN_in_out_L_out_tilde =  0.0, f_1PN_in_out_LL = 0.0, f_1PN_in_out_e_in = 0.0, f_1PN_in_out_i = 0.0;
    if (include_1PN_inner_outer_terms == TRUE)
    {
        f_1PN_in_out_m1_m2 = m1_times_m1 + m1_times_m2 + m2_times_m2;
        f_1PN_in_out_L_in_tilde = L_in/(m1_times_m2/m1_plus_m2);
        f_1PN_in_out_L_out_tilde = L_out/(m1_plus_m2*m3/m1_plus_m2_plus_m3);
        f_1PN_in_out_LL = f_1PN_in_out_L_in_tilde*f_1PN_in_out_L_out_tilde*m1_plus_m2*(4.0*m1_plus_m2 + 3.0*m3);
        f_1PN_in_out_e_in = (2.0 - 5.0*e_in_p2)*m1_times_m1 + 3.0*(-2.0 + e_in_p2)*m1_times_m1 + (2.0 - 5.0*e_in_p2)*m2_times_m2;
        f_1PN_in_out_i = 3.0*a_in*m1_plus_m2_plus_m3*cositot*(f_1PN_in_out_e_in + 3.0*e_in_p2*f_1PN_in_out_m1_m2*cos_2g_in);
    }
    	

//    printf("k_div_T %g %g %g \n",k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3);

	double R1_div_a_in = R1/a_in;
	double R1_div_a_in_p2 = R1_div_a_in*R1_div_a_in;
	double R1_div_a_in_p5 = R1_div_a_in_p2*R1_div_a_in_p2*R1_div_a_in;
	double R1_div_a_in_p6 = R1_div_a_in*R1_div_a_in_p5;
	double R2_div_a_in = R2/a_in;
	double R2_div_a_in_p2 = R2_div_a_in*R2_div_a_in;
	double R2_div_a_in_p5 = R2_div_a_in_p2*R2_div_a_in_p2*R2_div_a_in;
	double R2_div_a_in_p6 = R2_div_a_in*R2_div_a_in_p5;
	double R3_div_a_out = R3/a_out;
	double R3_div_a_out_p2 = R3_div_a_out*R3_div_a_out;
	double R3_div_a_out_p5 = R3_div_a_out_p2*R3_div_a_out_p2*R3_div_a_out;
	double R3_div_a_out_p6 = R3_div_a_out*R3_div_a_out_p5;
	
	double m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1 = m2_div_m1*R1_div_a_in_p6*k_div_T_tides_star1;
	double m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2 = m1_div_m2*R2_div_a_in_p6*k_div_T_tides_star2;
	double m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3 = m1_plus_m2_div_m3*R3_div_a_out_p6*k_div_T_tides_star3;
    
    double a_in_p2 = a_in*a_in;
    double a_in_p3 = a_in_p2*a_in;
    double a_in_p4 = a_in_p3*a_in;
	double n_in = sqrt(CONST_G*m1_plus_m2/a_in_p3); // mean orbital angular speed
    double spin_angular_frequency1_div_n_in = spin_angular_frequency1/n_in;
    double spin_angular_frequency2_div_n_in = spin_angular_frequency2/n_in;

    double a_out_p2 = a_out*a_out;
    double a_out_p3 = a_out_p2*a_out;
    double a_out_p4 = a_out_p3*a_out;
	double n_out = sqrt(CONST_G*m1_plus_m2_plus_m3/a_out_p3); // mean orbital angular speed
    double spin_angular_frequency3_div_n_out = spin_angular_frequency3/n_out;

    double f_tides1_in = f_tides1(e_in_p2);
    double f_tides2_in = f_tides2(e_in_p2);
    double f_tides3_in = f_tides3(e_in_p2);
    double f_tides4_in = f_tides4(e_in_p2);
    double f_tides5_in = f_tides5(e_in_p2);

    double f_tides1_out = f_tides1(e_out_p2);
    double f_tides2_out = f_tides2(e_out_p2);
    double f_tides3_out = f_tides3(e_out_p2);
    double f_tides4_out = f_tides4(e_out_p2);
    double f_tides5_out = f_tides5(e_out_p2);

    

    /************************************************
     * the calculations of the ODE right-hand-sides *
     * **********************************************/
     	

    /*******************************
     * e_in_dot                    *
     * *****************************/
    
	double e_in_dot_newtonian = 0.0;
    double e_in_dot_GR_1PN_in_out = 0.0;
    double e_in_dot_GR_25PN_in = 0.0;
    double e_in_dot_tides = 0.0;    

    /* Newtonian point particle -- up and including octupole order */    
    if (ignore_tertiary == FALSE)
    {
        e_in_dot_newtonian = C2*(l_in_p2/G_in)*(30.0*e_in*sinitot_p2*sin_2g_in) \
            + C3*e_out*(l_in_p2/G_in)*(35.0*cosphi*sinitot_p2*e_in_p2*sin_2g_in \
                - 10.0*cositot*sinitot_p2*cos_g_in*sin_g_out*l_in_p2 \
                - A*(sin_g_in*cos_g_out - cositot*cos_g_in*sin_g_out));
    }
    
    /* post-Newtonian point particle */
    if ((ignore_tertiary == FALSE) && (include_1PN_inner_outer_terms == TRUE))
    {
        e_in_dot_GR_1PN_in_out = -c_9div16*CONST_G*CONST_G*a_in*e_in*l_in*m1_times_m2*f_1PN_in_out_m1_m2*m3*sinitot_p2*sin_2g_in/( CONST_C_LIGHT_P2*a_out_p3*l_out_p3*L_in*m1_plus_m2*m1_plus_m2);
    }
    if (include_25PN_inner_terms == TRUE)
    {
        e_in_dot_GR_25PN_in = -c_304div15*CONST_G_P3*m1_times_m2*m1_plus_m2*(e_in/(CONST_C_LIGHT_P5*a_in_p4*l_in_p5))*f_25PN_e_in;
    }
    
    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        double e_in_dot_tides_star1 = -27.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency1_div_n_in);
        double e_in_dot_tides_star2 = -27.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency2_div_n_in);            

        if (e_in <= threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero)
        {
            e_in_dot_tides_star1 = e_in_dot_tides_star2 = 0.0;
//            printf("e_in_dot_tides %g\n",e_in_dot_tides_star1);
            
        }
        e_in_dot_tides = e_in_dot_tides_star1 + e_in_dot_tides_star2;
    }
    
    /* combined */

    double e_in_dot = 0.0;

    if (e_in <= threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero)
	{
//        printf("setting e_in_dot zero %g %g %g %g %g \n",e_in,e_in_dot_newtonian,e_in_dot_GR_1PN_in_out,e_in_dot_GR_25PN_in,e_in_dot_tides);
	    e_in_dot = 0.0;
	}
    else
	{
	    e_in_dot = e_in_dot_newtonian + e_in_dot_GR_1PN_in_out + e_in_dot_GR_25PN_in + e_in_dot_tides;
    }
	Ith(ydot,1) = -1.0*pow(10.0,-x)*e_in_dot/log(10.0);



    /*******************************
     * e_out_dot                       *
     * *****************************/

    double e_out_dot_newtonian = 0.0;
    double e_out_dot_GR_25PN_out = 0.0;
    double e_out_dot_tides = 0.0;
    
    if (ignore_tertiary == FALSE)
    {
        /* Newtonian point particle -- up and including octupole order */
        e_out_dot_newtonian = -C3*e_in*(l_out_p2/G_out)*(10.0*cositot*sinitot_p2*l_in_p2*sin_g_in*cos_g_out \
            + A*(cos_g_in*sin_g_out - cositot*sin_g_in*cos_g_out));
            
        /* post-Newtonian point particle */        
        if (include_25PN_outer_terms == TRUE)
        {
            e_out_dot_GR_25PN_out = -c_304div15*CONST_G_P3*m1_plus_m2*m3*m1_plus_m2_plus_m3*(e_out/(CONST_C_LIGHT_P5*a_out_p4*l_out_p5))*f_25PN_e_out;
        }
    
        /* tides */
        if (include_outer_tidal_terms == TRUE) // ad hoc/approximate: treats inner binary as point particle
        {
            double e_out_dot_tides_star3 = -27.0*(1.0+m1_plus_m2_div_m3)*m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3*R3_div_a_out_p2*e_out*pow(l_out,-13.0)*(f_tides3_out \
                - c_11div18*l_out_p3*f_tides4_out*spin_angular_frequency3_div_n_out);
            e_out_dot_tides = e_out_dot_tides_star3;
        }

        /* mass transfer */
        /* wind mass loss */
    }
    
    /* combined */
	double e_out_dot = e_out_dot_newtonian + e_out_dot_GR_25PN_out + e_out_dot_tides;
	Ith(ydot,2) = -1.0*pow(10.0,-y)*e_out_dot/log(10.0);



    /*******************************
     * g_in_dot                    *
     * *****************************/

	double g_in_dot_newtonian = 0.0;
    double g_in_dot_GR_1PN_in = 0.0;
    double g_in_dot_GR_1PN_in_out = 0.0;
    double g_in_dot_tides = 0.0;    
     
    /* Newtonian point particle -- up and including octupole order */
    if (ignore_tertiary == FALSE)
    {
        g_in_dot_newtonian = 6.0*C2*((1.0/G_in)*(4.0*cositot_p2 + (5.0*cos_2g_in - 1.0)*(l_in_p2 - cositot_p2)) \
                + (cositot/G_out)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in))) \
            - C3*e_out*(e_in*((1.0/G_out) + (cositot/G_in))*(sin_g_in*sin_g_out*(10.0*(3.0*cositot_p2 \
                - 1.0)*(1.0 - e_in_p2) + A) - 5.0*B*cositot*cosphi) \
                - (l_in_p2/(e_in*G_in))*(sin_g_in*sin_g_out*10.0*cositot*sinitot_p2*(1.0 - 3.0*e_in_p2) \
                + cosphi*(3.0*A - 10.0*cositot_p2 + 2.0)));
    }
    
    /* post-Newtonian point particle */  
    if (include_1PN_inner_terms == TRUE)
    {
        g_in_dot_GR_1PN_in = (3.0/(CONST_C_LIGHT_P2*a_in*l_in_p2))*pow(CONST_G*m1_plus_m2/a_in,3.0/2.0);
    }
    if ((ignore_tertiary == FALSE) && (include_1PN_inner_outer_terms == TRUE))
    {
        g_in_dot_GR_1PN_in_out = CONST_G*CONST_G*m1_times_m2_times_m3/(16.0*a_out_p3*CONST_C_LIGHT_P2*l_out_p3*m1_plus_m2*m1_plus_m2*m1_plus_m2_plus_m3)*( (a_in/G_in)*m1_plus_m2_plus_m3*( \
            l_in_p2*(5.0*m1_times_m2 - 3.0*m1_times_m2 + 5.0*m2_times_m2) - 9.0*f_1PN_in_out_m1_m2*( l_in_p2*cos_2g_in + 2.0*cositot_p2*sin_g_in*sin_g_in ) ) \
            + (1.0/G_out)*(-8.0*f_1PN_in_out_LL + f_1PN_in_out_i) );
    }
	
    /* tides/rotation */
    if (include_inner_tidal_terms == TRUE)
    {    
        double g_in_dot_tides_star1 = R1_div_a_in_p5*AMC_star1*(n_in/(l_in_p4))*(15.0*m2_div_m1*f_tides4_in/l_in_p6 \
		+ (1.0+m2_div_m1)*spin_angular_frequency1_div_n_in*spin_angular_frequency1_div_n_in);
        double g_in_dot_tides_star2 = R2_div_a_in_p5*AMC_star2*(n_in/(l_in_p4))*(15.0*m1_div_m2*f_tides4_in/l_in_p6 \
		+ (1.0+m1_div_m2)*spin_angular_frequency2_div_n_in*spin_angular_frequency2_div_n_in);
        g_in_dot_tides = g_in_dot_tides_star1 + g_in_dot_tides_star2;
    }
    
    /* combined */  
    double g_in_dot = g_in_dot_newtonian + g_in_dot_GR_1PN_in + g_in_dot_GR_1PN_in_out + g_in_dot_tides;
	Ith(ydot,3) = g_in_dot;


    /*******************************
     * g_out_dot                   *
     * *****************************/
     
    double g_out_dot_newtonian = 0.0;
    double g_out_dot_GR_1PN_out = 0.0;
    double g_out_dot_GR_1PN_in_out = 0.0;
    double g_out_dot_tides = 0.0;
    
    if (ignore_tertiary == FALSE)
    {
        /* Newtonian point particle -- up and including octupole order */
        g_out_dot_newtonian = 3.0*C2*((2.0*cositot/G_in)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)) \
                + (1.0/G_out)*(4.0 + 6.0*e_in_p2 + (5.0*cositot_p2 - 3.0)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)))) \
            + C3*e_in*(sin_g_in*sin_g_out*(((4.0*e_out_p2 + 1.0)/(e_out*G_out))*10.0*cositot*sinitot_p2*l_in_p2 \
                - e_out*((1.0/G_in) + (cositot/G_out))*(A + 10.0*(3.0*cositot_p2 - 1.0)*l_in_p2)) \
                + cosphi*(5.0*B*cositot*e_out*((1.0/G_in) + (cositot/G_out)) + ((4.0*e_out_p2 + 1.0)/(e_out*G_out))*A));
    
        /* post-Newtonian point particle */  
        if (include_1PN_outer_terms == TRUE)
        {    
            g_out_dot_GR_1PN_out = (3.0/(CONST_C_LIGHT_P2*a_out*l_out_p2))*pow(CONST_G*m1_plus_m2_plus_m3/a_out,3.0/2.0);
        }
        if (include_1PN_inner_outer_terms == TRUE)
        {
            g_out_dot_GR_1PN_in_out = CONST_G*CONST_G*m1_times_m2_times_m3/(16.0*a_out_p3*CONST_C_LIGHT_P2*l_out_p3*m1_plus_m2*m1_plus_m2*m1_plus_m2_plus_m3)*( (8.0*f_1PN_in_out_LL - f_1PN_in_out_i) \
                - (1.0/(2.0*G_out))*( 2.0*cositot*(-8.0*f_1PN_in_out_LL + f_1PN_in_out_i) - 16.0*m1_plus_m2*m1_plus_m2*cositot*f_1PN_in_out_L_in_tilde*f_1PN_in_out_L_out_tilde*( \
                    16.0*m1_plus_m2*f_1PN_in_out_L_in_tilde*f_1PN_in_out_L_out_tilde*(7.0*m1_plus_m2 + 6.0*m3)*cositot \
                        + c_3div2*a_in*m1_plus_m2_plus_m3*(-f_1PN_in_out_e_in*(1.0 + 3.0*(2.0*cositot_p2 - 1.0)) + 18.0*e_in_p2*f_1PN_in_out_m1_m2*cos_2g_in*sinitot_p2) ) ) );
        }
    
        /* tides */
        if (include_outer_tidal_terms == TRUE)
        {    
            double g_out_dot_tides_star3 = R3_div_a_out_p5*AMC_star3*(n_out/(l_out_p4))*(15.0*m1_plus_m2_div_m3*f_tides4_out/l_out_p6 \
            + (1.0+m1_plus_m2_div_m3)*spin_angular_frequency3_div_n_out*spin_angular_frequency3_div_n_out);
            g_out_dot_tides = g_out_dot_tides_star3;
        }
        
    }
    
    /* combined */
    double g_out_dot = g_out_dot_newtonian + g_out_dot_GR_1PN_out + g_out_dot_GR_1PN_in_out + g_out_dot_tides;
	Ith(ydot,4) = g_out_dot;
	
    
    
    /*******************************
     * h_in_dot                    *
     * *****************************/
     
    double h_in_dot_newtonian = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    /* 2013MNRAS.431.2155N Eq. B8 */
    if (ignore_tertiary == FALSE)
    {
        h_in_dot_newtonian = -3.0*C2*(G_tot/(G_in*G_out*sinitot))*(2.0 + 3.0*e_in_p2 - 5.0*e_in_p2*cos_2g_in)*2.0*sinitot*cositot \
            - C3*e_in*e_out*(G_tot/(G_in*G_out))*(5.0*B*cositot*cosphi - A*sin_g_in*sin_g_out + 10.0*(1.0 - 3.0*cositot_p2)*l_in_p2*sin_g_in*sin_g_out);
    }
    Ith(ydot,5) = h_in_dot_newtonian;
    Ith(ydot,6) = h_in_dot_newtonian;
   
   
   
    /*******************************
     * a_in_dot                    *
     * *****************************/

	double a_in_dot_GR_25PN_in = 0.0;
    double a_in_dot_tides = 0.0;    
    
    /* post-Newtonian point particle */  
    if (include_25PN_inner_terms == TRUE)
    {
        a_in_dot_GR_25PN_in = -c_64div5*CONST_G_P3*m1_times_m2*m1_plus_m2*(1.0/(CONST_C_LIGHT_P5*a_in_p3*l_in_p7))*f_25PN_a_in;
    }

    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        double a_in_dot_tides_star1 = -6.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency1_div_n_in);
        double a_in_dot_tides_star2 = -6.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency2_div_n_in);            
        a_in_dot_tides = a_in_dot_tides_star1 + a_in_dot_tides_star2;
        
        data->tidal_E1_dot = CONST_G*m1*m2/(2.0*a_in*a_in)*a_in_dot_tides_star1;
        data->tidal_E2_dot = CONST_G*m1*m2/(2.0*a_in*a_in)*a_in_dot_tides_star2;
    }

   
    /* combined */
    double a_in_dot = a_in_dot_GR_25PN_in + a_in_dot_tides;
	Ith(ydot,7) = a_in_dot;



    /*******************************
     * a_out_dot                   *
     * *****************************/

    double a_out_dot_GR_25PN_out = 0.0;
    double a_out_dot_tides = 0.0;
    
    if (ignore_tertiary == FALSE)
    {
        /* post-Newtonian point particle */  
        if (include_25PN_outer_terms == TRUE)
        {
            a_out_dot_GR_25PN_out = -c_64div5*CONST_G_P3*m1_plus_m2*m3*m1_plus_m2_plus_m3*(1.0/(CONST_C_LIGHT_P5*a_out_p3*l_out_p7))*f_25PN_a_out;
        }
    
        /* tides */
        if (include_outer_tidal_terms == TRUE)
        {
            double a_out_dot_tides_star3 = -6.0*(1.0+m1_plus_m2_div_m3)*m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3*R3_div_a_out_p2*a_out*pow(l_out,-15.0)*(f_tides1_out \
                - l_out_p3*f_tides2_out*spin_angular_frequency3_div_n_out);
            a_out_dot_tides = a_out_dot_tides_star3;
            
            data->tidal_E3_dot = CONST_G*(m1+m2)*m3/(2.0*a_out*a_out)*a_out_dot_tides_star3;
        }
        
    }
    
    /* combined */
    double a_out_dot = a_out_dot_GR_25PN_out + a_out_dot_tides;
	Ith(ydot,8) = a_out_dot;



    /**********************************************
     * cositot_dot                                *
     * due to dynamical triple interaction only!  *
     * ********************************************/

    double cositot_dot = 0.0;
    if (ignore_tertiary == FALSE)
    {
        double G_in_dot = -G_in*e_in*(e_in_dot_newtonian+e_in_dot_GR_1PN_in_out)/l_in_p2;
        double G_out_dot = -G_out*e_out*e_out_dot_newtonian/l_out_p2;
        cositot_dot = (-1.0/(G_in*G_out))*(G_in_dot*(G_in + G_out*cositot) + G_out_dot*(G_out + G_in*cositot));
    }
    
	Ith(ydot,9) = cositot_dot;



    /************************************
     * spin_angular_frequency1_dot      *
     * **********************************/

	double spin_angular_frequency1_dot_tides = 0.0;
    
    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        spin_angular_frequency1_dot_tides = 3.0*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*(m2_div_m1/(gyration_radius_star1*gyration_radius_star1)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency1_div_n_in);
    }

    double spin_angular_frequency1_dot = spin_angular_frequency1_dot_tides;
	Ith(ydot,10) = spin_angular_frequency1_dot;
    


    /************************************
     * spin_angular_frequency2_dot      *
     * **********************************/

	double spin_angular_frequency2_dot_tides = 0.0;
    
    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        spin_angular_frequency2_dot_tides = 3.0*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*(m1_div_m2/(gyration_radius_star2*gyration_radius_star2)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency2_div_n_in);
    }
    
    double spin_angular_frequency2_dot = spin_angular_frequency2_dot_tides;
	Ith(ydot,11) = spin_angular_frequency2_dot;


    /************************************
     * spin_angular_frequency3_dot      *
     * **********************************/

	double spin_angular_frequency3_dot_tides = 0.0;

    if (ignore_tertiary == FALSE)
    {
        /* tides */
        if (include_outer_tidal_terms == TRUE)
        {
            spin_angular_frequency3_dot_tides = 3.0*m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3*(m1_plus_m2_div_m3/(gyration_radius_star3*gyration_radius_star3)) \
                *(n_out/(l_out_p6*l_out_p6))*(f_tides2_out - l_out_p3*f_tides5_out*spin_angular_frequency3_div_n_out);
        }
    }
    
    double spin_angular_frequency3_dot = spin_angular_frequency3_dot_tides;
	Ith(ydot,12) = spin_angular_frequency3_dot;
    
	return 0;
}

/***************************************************
 * separate smaller functions for right-hand sides *
 ***************************************************/

/* tides (1981A&A....99..126H) */
double f_tides1(double e_p2)
{
    return 1.0 + e_p2*(c_31div2 + e_p2*(c_255div8 + e_p2*(c_185div16 + e_p2*c_25div64)));
}
double f_tides2(double e_p2)
{
    return 1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16));
}
double f_tides3(double e_p2)
{
    return 1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64));
}
double f_tides4(double e_p2)
{
    return 1.0 + e_p2*(c_3div2 + e_p2*c_1div8);
}
double f_tides5(double e_p2)
{
    return 1.0 + e_p2*(3.0 + e_p2*c_3div8);
}
double f_25PN_e(double e_p2)
{
    return 1.0 + e_p2*c_121div304;
}
double f_25PN_a(double e_p2)
{
    return 1.0 + e_p2*(c_73div24 + e_p2*c_37div96);
}


/* root finding functions */
int froot_delaunay(realtype t, N_Vector yev, realtype *gout, void *data_f)
{
	UserData data;
	data = (UserData) data_f;

    /* by default checks are not done (reduces computational cost) */
    /* setting gout[i] to >0 means no root will be found */
    gout[0] = 1e10;
    gout[1] = 1e10;
    gout[2] = 1e10;
    gout[3] = 1e10;
    gout[4] = 1e10;
    gout[5] = 1e10;

    bool check_for_dynamical_stability = data->check_for_dynamical_stability;
    bool check_for_inner_collision = data->check_for_inner_collision;
    bool check_for_outer_collision = data->check_for_outer_collision;
    bool check_for_inner_RLOF = data->check_for_inner_RLOF;
    bool check_for_outer_RLOF = data->check_for_outer_RLOF;

    double m1 = data->m1;
    double m2 = data->m2;				
    double m3 = data->m3;		
    
    double R1 = data->R1;
    double R2 = data->R2;
    double R3 = data->R3;

    double a_in = Ith(yev,7);
    double a_out = Ith(yev,8);
    double e_in = 1.0 - pow(10.0,Ith(yev,1)); /* Current e1 */
    double e_out = 1.0 - pow(10.0,Ith(yev,2)); /* Current e2 */

    double rp_in = a_in*(1.0 - e_in);
    double rp_out = a_out*(1.0 - e_out);
    
    if (check_for_dynamical_stability == TRUE)
    {
        /*	check for dynamical stability */
        /*  at the moment: uses Mardling & Aarseth criterion (2001MNRAS.321..398M) */
        /*  in future other criteria could be implemented as well */

        double itot = acos(Ith(yev,9));
        double a_out_div_a_in_crit = a_out_div_a_in_dynamical_stability(m1,m2,m3,e_out,itot);
        
        gout[0] = fabs(a_out/a_in - a_out_div_a_in_crit);
    }

    if (check_for_inner_collision == TRUE)
    {
        /*	check for collision at periastron (inner binary)	*/
        gout[1] = rp_in - (R1 + R2);
    }
    
    if (check_for_outer_collision == TRUE)
    {
        /*	check for "collision" at periastron (outer binary)	*/
        /*  the "radius" of the inner binary is set to the inner apocenter distance */
        double ra_in = a_in*(1.0 + e_in);
        gout[2] = rp_out - (R3 + ra_in);
    }
    
    if (check_for_inner_RLOF == TRUE)
    {
        double spin_angular_frequency1 = Ith(yev,10);
        double spin_angular_frequency2 = Ith(yev,11);
        double spin_angular_frequency_inner_orbit_periapse = sqrt( CONST_G*(m1+m2)*(1.0+e_in)/(rp_in*rp_in*rp_in) );
        double f1 = spin_angular_frequency1/spin_angular_frequency_inner_orbit_periapse;
        double f2 = spin_angular_frequency2/spin_angular_frequency_inner_orbit_periapse;
        
//        double roche_radius_pericenter_inner_star1 = roche_radius_pericenter_eggleton(rp_in, m1/m2);        
        double roche_radius_pericenter_inner_star1 = roche_radius_pericenter_sepinsky(rp_in, m1/m2, e_in, f1);
        gout[3] = R1 - roche_radius_pericenter_inner_star1;

//        double roche_radius_pericenter_inner_star2 = roche_radius_pericenter_eggleton(rp_in, m2/m1);        
        double roche_radius_pericenter_inner_star2 = roche_radius_pericenter_sepinsky(rp_in, m2/m1, e_in, f2);

        gout[4] = R2 - roche_radius_pericenter_inner_star2;
    }
    if (check_for_outer_RLOF == TRUE)
    {

        double spin_angular_frequency3 = Ith(yev,12);
        double spin_angular_frequency_outer_orbit_periapse = sqrt( CONST_G*(m1+m2+m3)*(1.0+e_out)/(rp_out*rp_out*rp_out) );        
        double f3 = spin_angular_frequency3/spin_angular_frequency_outer_orbit_periapse;
        
//        double roche_radius_pericenter_outer_star3 = roche_radius_pericenter_eggleton(rp_out, m3/(m1+m2));        
        double roche_radius_pericenter_outer_star3 = roche_radius_pericenter_sepinsky(rp_out, m3/(m1+m2), e_out, f3);
        gout[5] = R3 - roche_radius_pericenter_outer_star3;
    }            
    
	return 0;
}

double a_out_div_a_in_dynamical_stability(double m1, double m2, double m3, double e_out, double itot)
{
    /* wrapper used in interface.py */
    
    return a_out_div_a_in_dynamical_stability_mardling_aarseth_01(m1,m2,m3,e_out,itot);
}

double a_out_div_a_in_dynamical_stability_mardling_aarseth_01(double m1, double m2, double m3, double e_out, double itot)
{
    /* Mardling & Aarseth criterion (2001MNRAS.321..398M) 
     * including the `ad hoc' inclination factor
     * itot is assumed to be in radians! */
     
    double q_out = m3/(m1+m2);
    double a_out_div_a_in_crit = (1.0/(1.0-e_out))*2.8*pow( (1.0+q_out)*(1.0+e_out)/sqrt(1.0-e_out),c_2div5)*(1.0 - 0.3*itot/M_PI);
    
    return a_out_div_a_in_crit;
}

double roche_radius_pericenter_eggleton(double rp, double q)
{
    /* 2007ApJ...660.1624S Eqs. (45) */    
    /* q is defined as m_primary/m_secondary */
    double q_pow_one_third = pow(q,c_1div3);
    double q_pow_two_third = q_pow_one_third*q_pow_one_third;
    return rp*0.49*q_pow_two_third/(0.6*q_pow_two_third + log(1.0 + q_pow_one_third));
}
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f)
{
    /* 2007ApJ...660.1624S Eqs. (47)-(52) */
    double log_q = log10(q);
    double A = f*f*(1.0 + e); // assumes pericenter
    double log_A = log10(A);

    double R_L_pericenter_eggleton = roche_radius_pericenter_eggleton(rp,q);
    double ratio = 0.0; // this is R_L divided by R_L_pericenter_eggleton

    if (log_q < 0.0)
    {
        if (log_A <= -0.1)
        {
            double c = 0.5*(1.0+A) + log_q;
            ratio = 1.0 + 0.11*(1.0-A) - 0.05*(1.0-A)*exp(-c*c);
        }
        if ((log_A > -0.1) && (log_A < 0.2))
        {
            double g_0 = 0.9978 - 0.1229*log_A - 0.1273*log_A*log_A;
            double g_1 = 0.001 + 0.02556*log_A;
            double g_2 = 0.0004 + 0.0021*log_A;
            ratio = g_0 + g_1*log_q * g_2*log_q*log_q;
        }
        if (log_A >= 0.2)
        {
            double num_0 = 6.3014*pow(log_A,1.3643);
            double den_0 = exp(2.3644*pow(log_A,0.70748)) - 1.4413*exp(-0.0000184*pow(log_A,-4.5693));
            double i_0 = num_0/den_0;

            double den_1 = 0.0015*exp(8.84*pow(log_A,0.282)) + 15.78;
            double i_1 = log_A/den_1;

            double num_2 = 1.0 + 0.036*exp(8.01*pow(log_A,0.879));
            double den_2 = 0.105*exp(7.91*pow(log_A,0.879));
            double i_2 = num_2/den_2;

            double den_3 = 1.38*exp(-0.035*pow(log_A,0.76)) + 23.0*exp(-2.89*pow(log_A,0.76));
            double i_3 = 0.991/den_3;

            double c = log_q + i_3;
            ratio = i_0 + i_1*exp(-i_2*c*c);
        }
    }
    if (log_q >= 0.0)
    {
        if (log_A <= -0.1)
        {
            ratio = 1.226 - 0.21*A - 0.15*(1.0-A)*exp( (0.25*A - 0.3)*pow(log_q,1.55) );
        }
        if ((log_A > -0.1) && (log_A < 0.2))
        {
            double log_A_p2 = log_A*log_A;
            double h_0 = 1.0071 - 0.0907*log_A - 0.0495*log_A_p2;
            double h_1 = -0.004 - 0.163*log_A - 0.214*log_A_p2;
            double h_2 = 0.00022 - 0.0108*log_A - 0.02718*log_A_p2;
            ratio = h_0 + h_1*log_q + h_2*log_q*log_q;
        }
        if (log_A >= 0.2)
        {
            double num_0 = 1.895*pow(log_A,0.837);
            double den_0 = exp(1.636*pow(log_A,0.789)) - 1.0;
            double j_0 = num_0/den_0;

            double num_1 = 4.3*pow(log_A,0.98);
            double den_1 = exp(2.5*pow(log_A,0.66)) + 4.7;
            double j_1 = num_1/den_1;

            double den_2 = 8.8*exp(-2.95*pow(log_A,0.76)) + 1.64*exp(-0.03*pow(log_A,0.76));
            double j_2 = 1.0/den_2;

//            double j_3 = 0.256*exp(-1.33*pow(log_A,2.9))*( 5.5*exp(1.33*pow(log_A,2.9)) + 1.0 );
            double j_3 = 0.256*(5.5 + exp(-1.33*pow(log_A,2.9)));

            ratio = j_0 + j_1*exp(-j_2*pow(log_q,j_3));
            
//            printf("log_A %g\n",log_A);
//            printf("1 %g %g %g \n",num_0,den_0,j_0);
//            printf("2 %g %g %g \n",num_1,den_1,j_1);            
//            printf("2 %g %g %g \n",den_2,j_2,j_3);            
//            printf("ratio %g %g %g \n",ratio);            
        }
    }
    if (ratio == 0.0)
    {
        printf("unrecoverable error occurred in function roche_radius_pericenter_sepinsky in ODE_system.c\n");
        printf("log_q %g log_A %g ratio %g\n",log_q,log_A,ratio);
        printf("rp %g q %g e %g f %g\n",rp,q,e,f);
        exit(-1);
    }
    
    return ratio*R_L_pericenter_eggleton;
}
