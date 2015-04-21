/*	Worker code for SecularTriple, a secular triple gravitational dynamics code taking into account Newtonian, 1PN and 2.5PN terms	*/
/*	The relevant ODEs are solved consistently for each user supplied timestep using CVODE (Cohen & Hindmarsh 1996)	*/

#include "main_code.h"

//#include "cvode/cvode.h"				    /* prototypes for CVODE fcts., consts. */
//#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
//#include "cvode/cvode_dense.h"			    /* prototype for CVDense */
//#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
//#include "cvode/sundials_types.h"			/* definition of type realtype */

int fev_triad(realtype t, N_Vector yev, N_Vector ydot, void *data_f)
{
	
	UserData data;
	data = (UserData) data_f;
	
	/*	Constants which appear in the ODE right hand sides	*/
	double m1 = data->m1;					
	double m2 = data->m2;				
	double m3 = data->m3;				

    double f_quad = data->f_quad;
    double f_oct = data->f_oct;
    double f_tides = data->f_tides;
    double f_mass_transfer = data->f_mass_transfer;
    double f_1PN_in = data->f_1PN_in;
    double f_1PN_out = data->f_1PN_out;
    double f_1PN_in_out = data->f_1PN_in_out;
    double f_25PN_in = data->f_25PN_in;
    double f_25PN_out = data->f_25PN_out;
    
    double x_in_vec[3], x_out_vec[3];
    double e_in_vec[3], e_out_vec[3];    
    double h_in_vec[3], h_out_vec[3];
    int i=0;
    for (i=0; i<3; i++)
    {
        x_in_vec[i] = Ith(yev,i+1);
        x_out_vec[i] = Ith(yev,i+4);
        h_in_vec[i] = Ith(yev,i+7);
        h_out_vec[i] = Ith(yev,i+10);

        Ith(ydot,i+13) = 0.0; // reserved for stellar spin vectors
        Ith(ydot,i+16) = 0.0;
    }
    double x_in = norm3(x_in_vec);
    double x_out = norm3(x_out_vec);
    double e_in = 1.0 - exp(x_in);
    double e_out = 1.0 - exp(x_out);
    
    double e_in_vec_unit[3], e_out_vec_unit[3];
    double h_in_vec_unit[3], h_out_vec_unit[3];
    double h_in = norm3(h_in_vec);
    double h_out = norm3(h_out_vec);
    for (i=0; i<3; i++)
    {
        e_in_vec_unit[i] = x_in_vec[i]/x_in;
        e_out_vec_unit[i] = x_out_vec[i]/x_out;
        h_in_vec_unit[i] = h_in_vec[i]/h_in;
        h_out_vec_unit[i] = h_out_vec[i]/h_out;
    }
    
    double q_in_vec_unit[3], q_out_vec_unit[3];
    cross3(h_in_vec_unit,e_in_vec_unit,q_in_vec_unit);
    cross3(h_out_vec_unit,e_out_vec_unit,q_out_vec_unit);
    
    double e_in_p2 = e_in*e_in;
    double e_out_p2 = e_out*e_out;
    double l_in_p2 = 1.0 - e_in_p2;
    double l_out_p2 = 1.0 - e_out_p2;
    double l_in = sqrt(l_in_p2);
    double l_out = sqrt(l_out_p2);
    double l_in_p3 = l_in*l_in_p2;
    double l_out_p3 = l_out*l_out_p2;
    double a_in = h_in*h_in*(m1+m2)/( CONST_G*m1*m1*m2*m2*l_in_p2 );
    double a_out = h_out*h_out*(m1+m2+m3)/( CONST_G*(m1+m2)*(m1+m2)*m3*m3*l_out_p2 );
    
    double phi_0 = (3.0*CONST_G/4.0)*(a_in*a_in/(a_out*a_out*a_out))*(m1*m2*m3)/(m1+m2);
    double eps_oct = (25.0/16.0)*(a_in/a_out)*(e_out/(1.0-e_out_p2))*(m1-m2)/(m1+m2);
    double tau1 = h_in/(sqrt(1.0-e_in_p2)*phi_0);
    double tau_out = h_out/(sqrt(1.0-e_out_p2)*phi_0);
    
    /******************
     * tidal friction *
     * ****************/
    double tF1 = 0.0;
    double tF2 = 0.0;
    
    double e_in_p4 = e_in_p2*e_in_p2;
    double e_in_p6 = e_in_p2*e_in_p4;
    
    double V1 = (9.0/tF1)*0.0;
    double V2 = 0.0;
    double W1 = 0.0;
    double W2 = 0.0;
    double X1 = 0.0;
    double X2 = 0.0;
    double Y1 = 0.0;
    double Y2 = 0.0;
    double Z1 = 0.0;
    double Z2 = 0.0;

    /************
     * PN terms *
     * **********/
    double Z_1PN_in = (3.0*sqrt(CONST_G)*CONST_G)*(m1+m2)*sqrt(m1+m2)/( a_in*a_in*sqrt(a_in)*CONST_C_LIGHT_P2*l_in_p2 );
     
    /********************
     * right-hand sides *
     ********************/
    double de_in_vec_dt[3], de_out_vec_dt[3];
    double dh_in_vec_dt[3], dh_out_vec_dt[3];
    double ds1_vec_dt[3], ds2_vec_dt[3];
    
    double h_in_vec_unit_cross_h_out_vec_unit[3];
    double h_in_vec_unit_cross_e_out_vec_unit[3];
    double e_in_vec_unit_cross_h_out_vec_unit[3];
    double e_in_vec_unit_cross_e_out_vec_unit[3];
    
    // all dot products below involve unit vectors //
    double e_in_dot_h_out = dot3(e_in_vec_unit,h_out_vec_unit);
    double e_in_dot_e_out = dot3(e_in_vec_unit,e_out_vec_unit);
    double h_in_dot_h_out = dot3(h_in_vec_unit,h_out_vec_unit);
    double h_in_dot_e_out = dot3(h_in_vec_unit,e_out_vec_unit);
    
    double scalar_grad_e_in_phi_oct_h_in_cross_h_out = (-1.0/l_out_p3)*( 5.0*e_in_dot_h_out*e_in - 2.0*eps_oct*(7.0*e_in_dot_h_out*e_in*e_in_dot_e_out*e_in - l_in_p2*h_in_dot_h_out*h_in_dot_e_out) );
//    double scalar_grad_e1_phi_oct_h1_cross_e2 = (eps_oct/l2_p3)*( (c_1div5
    
    for (i=0; i<3; i++)
    {
        de_in_vec_dt[i] = (Z1 + Z2 + f_1PN_in*Z_1PN_in)*e_in*q_in_vec_unit[i] - (Y1 + Y2)*e_in*h_in_vec_unit[i] - (V1 + V2)*e_in_vec_unit[i] \
            + (1.0/tau1)*( l_in ) ;
    }
    
	/*	Some mass combinations	*/
//	double q_f = m2/m1;
//	double qinv_f = 1.0/q_f;
//	double m_prod_inner_f = m1_f*m2_f;
//	double m_tot_inner_f = m1_f+m2_f;
//	double m_prod_triple_f = m_prod_inner_f*m3_f;
//	double m_tot_triple_f = m_tot_inner_f+m3_f;
//	double m_prod_tot_inner_f = m_prod_inner_f*m_tot_inner_f;
//	double tildefme1_f = -10.0*m1_f*m1_f + 6.0*m1_f*m2_f - 10.0*m2_f*m2_f;

	/*	Eccentricity functions	*/
//	double e1_f2 = e1_f*e1_f;
//	double e1_f4 = e1_f2*e1_f2;
//	double e1_f2com = 1.0 - e1_f2;		/* 'com' stands for complement */

//	double e2_f2 = e2_f*e2_f;
//	double e2_f4 = e2_f2*e2_f2;
//	double e2_f2com = 1.0 - e2_f2;

//	double fme1_f = (2.0 - 5.0*e1_f2)*(m1_f*m1_f + m2_f*m2_f) - 3.0*(2.0 - e1_f2)*m1_f*m2_f;
//	double f_GR_adot1 = 1.0 + (73.0/24.0)*e1_f2 + (37.0/96.0)*e1_f4;
//	double f_GR_edot1 = 1.0 + (121.0/304.0)*e1_f2;

//	double f_GR_adot2 = 1.0 + (73.0/24.0)*e2_f2 + (37.0/96.0)*e2_f4;
//	double f_GR_edot2 = 1.0 + (121.0/304.0)*e2_f2;
}	

int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data_f)
{

    /* TODO: optimalizations: pre-load fractions */
    
	UserData data;
	data = (UserData) data_f;
	
	/*	constants which appear in the ODE right hand sides	*/
	double m1 = data->m1; // mass
	double m2 = data->m2;				
	double m3 = data->m3;		

    double R1 = data->R1; // radius
    double R2 = data->R2;
    double R3 = data->R3;

    double AMC_star1 = data->AMC_star1; // Apsidal Motion Constant
    double AMC_star2 = data->AMC_star2; // Apsidal Motion Constant
    double AMC_star3 = data->AMC_star3; // Apsidal Motion Constant  
    
    double gyration_radius_star1 = data->gyration_radius_star1; // gyration radius (NOT squared)     
    double gyration_radius_star2 = data->gyration_radius_star2; // gyration radius (NOT squared)     
    double gyration_radius_star3 = data->gyration_radius_star3; // gyration radius (NOT squared)             

    double k_div_T_tides_star1 = data->k_div_T_tides_star1; // tidal dissipation constant
    double k_div_T_tides_star2 = data->k_div_T_tides_star2;    
    double k_div_T_tides_star3 = data->k_div_T_tides_star3;

    double m1dot = data->m1dot; // mass transfer rate
    double gamma_in = data->gamma_in; // 
    double gamma_out = data->gamma_out;
    double mu_in = data->mu_in; //
    double mu_out = data->mu_out;
    
    double f_quad = data->f_quad;
    double f_oct = data->f_oct;
    double f_tides = data->f_tides;
    double f_mass_transfer = data->f_mass_transfer;
    double f_1PN_in = data->f_1PN_in;
    double f_1PN_out = data->f_1PN_out;
    double f_1PN_in_out = data->f_1PN_in_out;
    double f_25PN_in = data->f_25PN_in;
    double f_25PN_out = data->f_25PN_out;

	/*	the ODE variables	*/
	double x = Ith(yev,1);
	double y = Ith(yev,2);
	double g_in = Ith(yev,3);
	double g_out = Ith(yev,4);
	double a_in = Ith(yev,5);
	double a_out = Ith(yev,6);
	double cositot = Ith(yev,7);
	double spin_angular_frequency1 = Ith(yev,8);
	double spin_angular_frequency2 = Ith(yev,9);

	double e_in = 1.0 - pow(10.0,x);
	double e_out = 1.0 - pow(10.0,y);

	/*	some mass combinations	*/
	double m1_div_m2 = m1/m2;
	double m2_div_m1 = 1.0/m1_div_m2;
    
    /* below: should be renamed/cleaned up */
	double m_prod_inner = m1*m2;
	double m_tot_inner = m1+m2;
	double m_prod_triple = m_prod_inner*m3;
	double m_tot_triple = m_tot_inner+m3;
	double m_prod_tot_inner = m_prod_inner*m_tot_inner;
	double tildefme1 = -10.0*m1*m1 + 6.0*m1*m2 - 10.0*m2*m2;

    double m1_plus_m2 = m1+m2;
    double m1_plus_m2_plus_m3 = m1_plus_m2+m3;
    double m1_times_m2 = m1*m2;
    double m1_times_m2_times_m3 = m1_times_m2*m3;
    
	/*	eccentricity functions	*/
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

//	double f_GR_adot1 = 1.0 + (73.0/24.0)*e1_p2 + (37.0/96.0)*e1_p4;
//	double f_GR_edot1 = 1.0 + (121.0/304.0)*e1_p2;

//	double f_GR_adot2 = 1.0 + (73.0/24.0)*e2_p2 + (37.0/96.0)*e2_p4;
//	double f_GR_edot2 = 1.0 + (121.0/304.0)*e2_p2;
	
    /* tides quantities */
	double R1_div_a_in = R1/a_in;
	double R1_div_a_in_p2 = R1_div_a_in*R1_div_a_in;
	double R1_div_a_in_p5 = R1_div_a_in_p2*R1_div_a_in_p2*R1_div_a_in;
	double R1_div_a_in_p6 = R1_div_a_in*R1_div_a_in_p5;
	double R2_div_a_in = R2/a_in;
	double R2_div_a_in_p2 = R2_div_a_in*R2_div_a_in;
	double R2_div_a_in_p5 = R2_div_a_in_p2*R2_div_a_in_p2*R2_div_a_in;
	double R2_div_a_in_p6 = R2_div_a_in*R2_div_a_in_p5;
	
	double m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1 = m2_div_m1*R1_div_a_in_p6*k_div_T_tides_star1;
	double m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2 = m1_div_m2*R2_div_a_in_p6*k_div_T_tides_star2;
    
    double a_in_p2 = a_in*a_in;
    double a_in_p3 = a_in_p2*a_in;
    double a_in_p4 = a_in_p3*a_in;
	double n_in = sqrt(CONST_G*m1_plus_m2/a_in_p3); // mean orbital angular speed
    double spin_angular_frequency1_div_n_in = spin_angular_frequency1/n_in;
    double spin_angular_frequency2_div_n_in = spin_angular_frequency2/n_in;

    double a_out_p2 = a_out*a_out;
    double a_out_p3 = a_out_p2*a_out;
    double a_out_p4 = a_out_p3*a_out;

    double f_tides1_in = f_tides1(e_in_p2);
    double f_tides2_in = f_tides2(e_in_p2);
    double f_tides3_in = f_tides3(e_in_p2);
    double f_tides4_in = f_tides4(e_in_p2);
    double f_tides5_in = f_tides5(e_in_p2);
                    
    /* PN quantities */
    double f_25PN_e_in = f_25PN_e(e_in_p2);
    double f_25PN_a_in = f_25PN_a(e_in_p2);

    double f_25PN_e_out = f_25PN_e(e_out_p2);
    double f_25PN_a_out = f_25PN_a(e_out_p2);
    double fme1 = (2.0 - 5.0*e_in_p2)*(m1*m1 + m2*m2) - 3.0*(2.0 - e_in_p2)*m1*m2;

    /* mass transfer quantities */
    double r_A1_P_div_a_in,r_A2_P_div_a_in,cos_Phi_P,X_L1_P_approx;
    if (f_mass_transfer != 0.0)
    {
        double f=0.0;
        X_L1_P_approx = 0.529 + 0.231*log(m1_div_m2) - f*f*(0.031 + 0.025*e_in)*(1.0 + 0.4*log(m1_div_m2)); // 2007ApJ...667.1170S Eq. A15
        r_A1_P_div_a_in = X_L1_P_approx*(1.0-e_in);
        
        r_A2_P_div_a_in = R2/a_in;
        r_A2_P_div_a_in = 0.0;
        cos_Phi_P = 1.0;        
    }

	/*	triple quantities */
    /* 2000ApJ...535..385F */
	double L_in = m1_times_m2*sqrt(CONST_G*a_in/(m1_plus_m2));
	double L_out = (m1_plus_m2)*m3*sqrt(CONST_G*a_out/(m1_plus_m2_plus_m3));
	double G_in = L_in*sqrt(l_in_p2);
	double G_out = L_out*sqrt(l_out_p2);

	double a_in_div_a_out = a_in/a_out;
	double C2 = CONST_G*c_1div16*(m1_times_m2_times_m3/m1_plus_m2)*pow(l_out,-3.0)*a_in_div_a_out*a_in_div_a_out/a_out;
	double C3 = -CONST_G*c_15div16*c_1div4*(m1_times_m2_times_m3/(m1_plus_m2*m1_plus_m2))*(m1-m2)*pow(l_out,-5.0)*a_in_div_a_out*a_in_div_a_out*a_in_div_a_out/a_out;
    C2 *= f_quad;
    C3 *= f_oct;

	if (cositot > 1.0)
	{
		cositot = 2.0 - cositot;
	}

	if (cositot < -1.0)
	{
		cositot = -2.0 - cositot;
	}

	double cositot_p2 = cositot*cositot;
	double sinitot = sqrt(1.0 - cositot_p2);
	double sinitot_p2 = sinitot*sinitot;

	double sin_g_in = sin(g_in);
	double sin_2g_in = sin(2.0*g_in);	/* 'd' stands for double angle */
	double sin_g_out = sin(g_out);
	double cos_g_in = cos(g_in);
	double cos_2g_in = cos(2.0*g_in);
	double cos_g_out = cos(g_out);
	
	/*	Required for octupole order terms	*/
	double B = 2.0 + 5.0*e_in_p2 - 7.0*e_in_p2*cos_2g_in;
	double A = 4.0 + 3.0*e_in_p2 - c_5div2*B*sinitot_p2;
	double cosphi = -cos_g_in*cos_g_out - cositot*sin_g_in*sin_g_out;
	

	/* ================================================================
	* The calculations of the ODE right hand sides
	* References: Ford et al. (2000); Blaes et al. (2002); Naoz et al. (2012)
	*/
	

    /*******************************
     * e_in_dot                       *
     * *****************************/
    
    /* Newtonian point particle -- up and including octupole order */
	double e_in_dot_newtonian = C2*(l_in_p2/G_in)*(30.0*e_in*sinitot_p2*sin_2g_in) \
		+ C3*e_out*(l_in_p2/G_in)*(35.0*cosphi*sinitot_p2*e_in_p2*sin_2g_in \
			- 10.0*cositot*sinitot_p2*cos_g_in*sin_g_out*l_in_p2 \
			- A*(sin_g_in*cos_g_out - cositot*cos_g_in*sin_g_out));
    
    /* post-Newtonian point particle */
    double e_in_dot_GR_1PN_in_out=0.0,e_in_dot_GR_25PN_in=0.0;
    if (f_1PN_in_out != 0.0)
    {
        e_in_dot_GR_1PN_in_out = 0.0;
        // -c_9div16*(CONST_G*CONST_G*pow(CONST_C_LIGHT,-2.0))*a_in*e_in*sqrt(l_in_p2)*m1*m2*(m1*m1 + m1*m2 + m2*m2) \
            *m3*sinitot_p2*sin_2g_in*pow( pow(a_out,3.0)*pow(l_out,3.0)*L_in*pow(m1+m2,2.0),-1.0); //
            // need to check this 
    }
    if (f_25PN_in != 0.0)
    {
        e_in_dot_GR_25PN_in = -c_304div15*CONST_G_P3*m1_times_m2*m1_plus_m2*(e_in/(CONST_C_LIGHT_P5*a_in_p4*l_in_p5))*f_25PN_e_in;
    }
    
    /* tides */
    double e_in_dot_tides_star1=0,e_in_dot_tides_star2=0;
    if (f_tides != 0.0)
    {
        e_in_dot_tides_star1 = -27.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency1_div_n_in);
        e_in_dot_tides_star2 = -27.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency2_div_n_in);            
    }
    
    /* mass transfer */
    double e_in_dot_mass_transfer=0.0;
    if (f_mass_transfer != 0.0)
    {
        e_in_dot_mass_transfer = (l_in/(2.0*M_PI))*(m1dot/m1)*( gamma_in*m1_div_m2*r_A2_P_div_a_in*cos_Phi_P + r_A1_P_div_a_in \
            + 2.0*(gamma_in*m1_div_m2-1.0)*(1.0-e_in) + 2.0*(1.0-gamma_in)*(mu_in + c_1div2)*(1.0-e_in*m1_div_m2/(1.0+m1_div_m2) ) );
//        if (e1 <= 0.0)
//        {
//            e1dot_mass_transfer = 0.0;
//        }
    }

    /* combined */
	double e_in_dot = e_in_dot_newtonian + f_1PN_in_out*e_in_dot_GR_1PN_in_out + f_25PN_in*e_in_dot_GR_25PN_in + f_tides*(e_in_dot_tides_star1 + e_in_dot_tides_star2) \
        + f_mass_transfer*e_in_dot_mass_transfer;
	Ith(ydot,1) = -1.0*pow(10.0,-x)*e_in_dot/log(10.0);


    /*******************************
     * e_out_dot                       *
     * *****************************/

    /* Newtonian point particle -- up and including octupole order */
	double e_out_dot_newtonian = -C3*e_in*(l_out_p2/G_out)*(10.0*cositot*sinitot_p2*l_in_p2*sin_g_in*cos_g_out \
		+ A*(cos_g_in*sin_g_out - cositot*sin_g_in*cos_g_out));
        
    /* post-Newtonian point particle */        
	double e_out_dot_GR_25PN_out=0.0;
    if (f_25PN_out != 0.0)
    {
        e_out_dot_GR_25PN_out = -c_304div15*CONST_G_P3*m1_plus_m2*m3*m1_plus_m2_plus_m3*(e_out/(CONST_C_LIGHT_P5*a_out_p4*l_out_p5))*f_25PN_e_out;
    }

    /* combined */
	double e_out_dot = e_out_dot_newtonian + f_25PN_out*e_out_dot_GR_25PN_out;
	Ith(ydot,2) = -1.0*pow(10.0,-y)*e_out_dot/log(10.0);


    /*******************************
     * g_in_dot                    *
     * *****************************/
     
    /* Newtonian point particle -- up and including octupole order */
	double g_in_dot_newtonian = 6.0*C2*((1.0/G_in)*(4.0*cositot_p2 + (5.0*cos_2g_in - 1.0)*(l_in_p2 - cositot_p2)) \
			+ (cositot/G_out)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in))) \
		- C3*e_out*(e_in*((1.0/G_out) + (cositot/G_in))*(sin_g_in*sin_g_out*(10.0*(3.0*cositot_p2 \
			- 1.0)*(1.0 - e_in_p2) + A) - 5.0*B*cositot*cosphi) \
			- (l_in_p2/(e_in*G_in))*(sin_g_in*sin_g_out*10.0*cositot*sinitot_p2*(1.0 - 3.0*e_in_p2) \
			+ cosphi*(3.0*A - 10.0*cositot_p2 + 2.0)));

    /* post-Newtonian point particle */  
    double g_in_dot_GR_1PN_in=0.0,g_in_dot_GR_1PN_in_out=0.0;
    if (f_1PN_in != 0.0)
    {
        g_in_dot_GR_1PN_in = (3.0/(CONST_C_LIGHT_P2*a_in*l_in_p2))*pow(CONST_G*m1_plus_m2/a_in,3.0/2.0);
    }
    if (f_1PN_in_out != 0.0)
    {
        g_in_dot_GR_1PN_in_out = 0.0;
        // = (CONST_G*pow(CONST_C_LIGHT,-2.0))*pow( 4.0*pow(a_out,3.0)*pow(l_out_p2,3.0/2.0)*(m1+m2),-1.0) \
            *( G_out*( 8.0*(m1+m2) + 6.0*m3)*cositot \
                - (l_in_p2/G_in)*(CONST_G)*a_in*m1*m2*m3*pow(8.0*(m1+m2),-1.0)*(tildefme1*(1.0 - 3.0*cositot_p2) \
                    + 18.0*(m1*m1 + m1*m2 + m2*m2)*cos_2g_in*sinitot_p2) \
                - ( (cositot/G1) + (1.0/G2) )*( G1*G2*(8.0*(m1+m2) + 6.0*m3) \
                    + cositot*a1*CONST_G*m1*m2*m3*pow(8.0*(m1+m2),-1.0) \
                    *(6.0*fme1 + 18.0*e1_p2*(m1*m1 + m1*m2 + m2*m2)*cosg1d) ) ); //
                    // need to check this
    }
	
    /* tides */
    double g_in_dot_tides_star1=0,g_in_dot_tides_star2=0;
    if (f_tides != 0.0)
    {    
        g_in_dot_tides_star1 = R1_div_a_in_p5*AMC_star1*(n_in/(l_in_p4))*(15.0*m2_div_m1*f_tides4_in/l_in_p6 \
		+ (1.0+m2_div_m1)*spin_angular_frequency1_div_n_in*spin_angular_frequency1_div_n_in);
        g_in_dot_tides_star2 = R2_div_a_in_p5*AMC_star2*(n_in/(l_in_p4))*(15.0*m1_div_m2*f_tides4_in/l_in_p6 \
		+ (1.0+m1_div_m2)*spin_angular_frequency2_div_n_in*spin_angular_frequency2_div_n_in);
    }
    /* mass transfer */    
    double g_in_dot_mass_transfer=0.0;
    if (f_mass_transfer != 0.0)
    {
        g_in_dot_mass_transfer = 0.0;
    }
    
    /* combined */    
	Ith(ydot,3) = g_in_dot_newtonian + f_1PN_in*g_in_dot_GR_1PN_in + f_1PN_in_out*g_in_dot_GR_1PN_in_out + f_tides*(g_in_dot_tides_star1 + g_in_dot_tides_star2) \
        + f_mass_transfer*g_in_dot_mass_transfer;


    /*******************************
     * g_out_dot                   *
     * *****************************/
     
    /* Newtonian point particle -- up and including octupole order */
	double g_out_dot_newtonian = 3.0*C2*((2.0*cositot/G_in)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)) \
			+ (1.0/G_out)*(4.0 + 6.0*e_in_p2 + (5.0*cositot_p2 - 3.0)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)))) \
		+ C3*e_in*(sin_g_in*sin_g_out*(((4.0*e_out_p2 + 1.0)/(e_out*G_out))*10.0*cositot*sinitot_p2*l_in_p2 \
			- e_out*((1.0/G_in) + (cositot/G_out))*(A + 10.0*(3.0*cositot_p2 - 1.0)*l_in_p2)) \
			+ cosphi*(5.0*B*cositot*e_out*((1.0/G_in) + (cositot/G_out)) + ((4.0*e_out_p2 + 1.0)/(e_out*G_out))*A));

    /* post-Newtonian point particle */  
    double g_out_dot_GR_1PN_out=0.0,g_out_dot_GR_1PN_in_out=0.0;
    if (f_1PN_in != 0.0)
    {    
        g_out_dot_GR_1PN_out = (3.0/(CONST_C_LIGHT_P2*a_out*l_out_p2))*pow(CONST_G*m1_plus_m2_plus_m3/a_out,3.0/2.0);
    }
    if (f_1PN_in_out != 0.0)
    {
        g_out_dot_GR_1PN_in_out = 0.0;
        //(CONST_G*pow(CONST_C_LIGHT_P2,-1.0))*pow( 4.0*pow(a2,3.0)*pow(l2_p2,3.0/2.0)*(m1+m2),-1.0) \
            *( -2.0*G1*(8.0*(m1+m2)+6.0*m3)*cositot - 3.0*CONST_G*a1*m1*m2*m3*pow(8.0*(m1+m2)*G2,-1.0) \
                    *( fme1 - 3.0*fme1*cositot_p2 + 9.0*e1_p2*(m1*m1+m1*m2+m2*m2)*cosg1d*sinitot_p2 ) \
                + ( (cositot/G2) + (1.0/G1) )*(G1*G2*(8.0*(m1+m2)+6.0*m3) \
                    + CONST_G*a1*m1*m2*m3*pow(8.0*(m1+m2),-1.0)*cositot*(6.0*fme1 \
                        + 18.0*e1_p2*(m1*m1+m1*m2+m2*m2)*cosg1d)) ); //
                        // need to check this
    }
    
    /* combined */
	Ith(ydot,4) = g_out_dot_newtonian + f_1PN_out*g_out_dot_GR_1PN_out + f_1PN_in_out*g_out_dot_GR_1PN_in_out;
	
   
    /*******************************
     * a_in_dot                    *
     * *****************************/

    /* post-Newtonian point particle */  
	double a_in_dot_GR_25PN_in=0.0;
    if (f_25PN_in != 0.0)
    {
        a_in_dot_GR_25PN_in = -c_64div5*CONST_G_P3*m1_times_m2*m1_plus_m2*(1.0/(CONST_C_LIGHT_P5*a_in_p3*l_in_p7))*f_25PN_a_in;
    }

    /* tides */
    double a_in_dot_tides_star1=0,a_in_dot_tides_star2=0;
    if (f_tides != 0.0)
    {
        a_in_dot_tides_star1 = -6.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency1_div_n_in);
        a_in_dot_tides_star2 = -6.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency2_div_n_in);            
    }

    /* mass transfer */
    double a_in_dot_mass_transfer=0.0;
    if (f_mass_transfer != 0.0)
    {
        a_in_dot_mass_transfer = (a_in/M_PI)*(m1dot/m1)*(1.0/l_in)*( e_in*r_A1_P_div_a_in + gamma_in*m1_div_m2*e_in*r_A2_P_div_a_in*cos_Phi_P \
            + (gamma_in*m1_div_m2-1.0)*(1.0 - e_in_p2) \
            + (1.0-gamma_in)*(mu_in + c_1div2)*l_in_p2*m1_div_m2/(1.0 + m1_div_m2) );
    }
        
    /* combined */
	Ith(ydot,5) = f_25PN_in*a_in_dot_GR_25PN_in + f_tides*(a_in_dot_tides_star1 + a_in_dot_tides_star2) + f_mass_transfer*a_in_dot_mass_transfer;


    /*******************************
     * a_out_dot                   *
     * *****************************/

    /* post-Newtonian point particle */  
	double a_out_dot_GR_25PN_out=0.0;
    if (f_25PN_out != 0.0)
    {
        a_out_dot_GR_25PN_out = -c_64div5*CONST_G_P3*m1_plus_m2*m3*m1_plus_m2_plus_m3*(1.0/(CONST_C_LIGHT_P5*a_out_p3*l_out_p7))*f_25PN_a_out;
    }

    /* combined */
	Ith(ydot,6) = f_25PN_out*a_out_dot_GR_25PN_out;


    /************************************
     * cositot_ot                       *
     * due to triple interaction only!  *
     * **********************************/

	double G_in_dot = (-1.0)*G_in*e_in*(e_in_dot_newtonian+f_1PN_in_out*e_in_dot_GR_1PN_in_out)/l_in_p2;
	double G_out_dot = (-1.0)*G_out*e_out*e_out_dot_newtonian/l_out_p2;
	double cositot_dot = (-1.0/(G_in*G_out))*(G_in_dot*(G_in + G_out*cositot) + G_out_dot*(G_out + G_in*cositot));

	Ith(ydot,7) = cositot_dot;


    /************************************
     * spin_angular_frequency1_dot      *
     * **********************************/

	double spin_angular_frequency1_dot=0.0;
    if (f_tides != 0.0)
    {
        spin_angular_frequency1_dot = 3.0*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*(m2_div_m1/(gyration_radius_star1*gyration_radius_star1)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency1_div_n_in);
    }
	Ith(ydot,8) = f_tides*spin_angular_frequency1_dot;
    

    /************************************
     * spin_angular_frequency2_dot      *
     * **********************************/

	double spin_angular_frequency2_dot=0.0;
    if (f_tides != 0.0)
    {
        spin_angular_frequency2_dot = 3.0*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*(m1_div_m2/(gyration_radius_star2*gyration_radius_star2)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency2_div_n_in);
    }
	Ith(ydot,9) = f_tides*spin_angular_frequency2_dot;

	return 0;
}

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

int froot(realtype t, N_Vector yev, realtype *gout, void *data)
{
	UserData data_f;
	data_f = (UserData) data;

/*	Check for dynamical stability (Marding & Aarseth 2001) */
	double e1_t = 1.0 - pow(10.0,Ith(yev,1)); /* Current e1 */
	double e2_t = 1.0 - pow(10.0,Ith(yev,2)); /* Current e2 */
	double q_out_t = data_f->m3/(data_f->m1 + data_f->m2);
	double itot_t = acos(Ith(yev,7));
	double beta_crit = (1.0/(1.0-e2_t))*2.8*pow( (1.0+q_out_t)*(1.0+e2_t)/sqrt(1.0-e2_t),2.0/5.0)*(1.0 \
		- 0.3*itot_t/M_PI);

	gout[0] = fabs(Ith(yev,6)/Ith(yev,5) - beta_crit);

/*	Check for collision at periastron (inner binary)	*/
	double a1_t = Ith(yev,5);
	double R1_t = data_f->R1;
	double R2_t = data_f->R2;
	gout[1] = a1_t*(1.0-e1_t) - (R1_t + R2_t);

	return 0;
}
