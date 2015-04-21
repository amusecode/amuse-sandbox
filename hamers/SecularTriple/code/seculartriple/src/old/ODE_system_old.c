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

    double f_oct = data->f_oct;
    double f_1PN_1 = data->f_1PN_1;
    double f_1PN_2 = data->f_1PN_2;
    double f_1PN_12 = data->f_1PN_12;
    double f_25PN_1 = data->f_25PN_1;
    double f_25PN_2 = data->f_25PN_2;
    
    double x1_vec[3], x2_vec[3];
    double e1_vec[3], e2_vec[3];    
    double h1_vec[3], h2_vec[3];
    int i=0;
    for (i=0; i<3; i++)
    {
        x1_vec[i] = Ith(yev,i+1);
        x2_vec[i] = Ith(yev,i+4);
        h1_vec[i] = Ith(yev,i+7);
        h2_vec[i] = Ith(yev,i+10);

        Ith(ydot,i+13) = 0.0; // reserved for stellar spin vectors
        Ith(ydot,i+16) = 0.0;
    }
    double x1 = norm3(x1_vec);
    double x2 = norm3(x2_vec);
    double e1 = 1.0 - exp(x1);
    double e2 = 1.0 - exp(x2);
    
    double e1_vec_unit[3], e2_vec_unit[3];
    double h1_vec_unit[3], h2_vec_unit[3];
    double h1 = norm3(h1_vec);
    double h2 = norm3(h2_vec);
    for (i=0; i<3; i++)
    {
        e1_vec_unit[i] = x1_vec[i]/x1;
        e2_vec_unit[i] = x2_vec[i]/x2;
        h1_vec_unit[i] = h1_vec[i]/h1;
        h2_vec_unit[i] = h2_vec[i]/h2;
    }
    
    double q1_vec_unit[3], q2_vec_unit[3];
    cross3(h1_vec_unit,e1_vec_unit,q1_vec_unit);
    cross3(h2_vec_unit,e2_vec_unit,q2_vec_unit);
    
    double e1_p2 = e1*e1;
    double e2_p2 = e2*e2;
    double l1_p2 = 1.0 - e1_p2;
    double l2_p2 = 1.0 - e2_p2;
    double l1 = sqrt(l1_p2);
    double l2 = sqrt(l2_p2);
    double l1_p3 = l1*l1_p2;
    double l2_p3 = l2*l2_p2;
    double a1 = h1*h1*(m1+m2)/( CONST_G*m1*m1*m2*m2*l1_p2 );
    double a2 = h2*h2*(m1+m2+m3)/( CONST_G*(m1+m2)*(m1+m2)*m3*m3*l2_p2 );
    
    double phi_0 = (3.0*CONST_G/4.0)*(a1*a1/(a2*a2*a2))*(m1*m2*m3)/(m1+m2);
    double eps_oct = (25.0/16.0)*(a1/a2)*(e2/(1.0-e2_p2))*(m1-m2)/(m1+m2);
    double tau1 = h1/(sqrt(1.0-e1_p2)*phi_0);
    double tau2 = h2/(sqrt(1.0-e2_p2)*phi_0);
    
    /******************
     * tidal friction *
     * ****************/
    double tF1 = 0.0;
    double tF2 = 0.0;
    
    double e1_p4 = e1_p2*e1_p2;
    double e1_p6 = e1_p2*e1_p4;
    
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
    double Z_1PN_1 = (3.0*sqrt(CONST_G)*CONST_G)*(m1+m2)*sqrt(m1+m2)/( a1*a1*sqrt(a1)*CONST_C_LIGHT_P2*l1_p2 );
     
    /********************
     * right-hand sides *
     ********************/
    double de1_vec_dt[3], de2_vec_dt[3];
    double dh1_vec_dt[3], dh2_vec_dt[3];
    double ds1_vec_dt[3], ds2_vec_dt[3];
    
    double h1_vec_unit_cross_h2_vec_unit[3];
    double h1_vec_unit_cross_e2_vec_unit[3];
    double e1_vec_unit_cross_h2_vec_unit[3];
    double e1_vec_unit_cross_e2_vec_unit[3];
    
    // all dot products below involve unit vectors //
    double e1_dot_h2 = dot3(e1_vec_unit,h2_vec_unit);
    double e1_dot_e2 = dot3(e1_vec_unit,e2_vec_unit);
    double h1_dot_h2 = dot3(h1_vec_unit,h2_vec_unit);
    double h1_dot_e2 = dot3(h1_vec_unit,e2_vec_unit);
    
    double scalar_grad_e1_phi_oct_h1_cross_h2 = (-1.0/l2_p3)*( 5.0*e1_dot_h2*e1 - 2.0*eps_oct*(7.0*e1_dot_h2*e1*e1_dot_e2*e1 - l1_p2*h1_dot_h2*h1_dot_e2) );
//    double scalar_grad_e1_phi_oct_h1_cross_e2 = (eps_oct/l2_p3)*( (c_1div5
    
    for (i=0; i<3; i++)
    {
        de1_vec_dt[i] = (Z1 + Z2 + f_1PN_1*Z_1PN_1)*e1*q1_vec_unit[i] - (Y1 + Y2)*e1*h1_vec_unit[i] - (V1 + V2)*e1_vec_unit[i] \
            + (1.0/tau1)*( l1 ) ;
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
	double m1 = data->m1;					
	double m2 = data->m2;				
	double m3 = data->m3;		

    double R1 = data->R1;
    double R2 = data->R2;
    double R3 = data->R3;

    double k_div_T_tides_star1 = data->k_div_T_tides_star1;
    double k_div_T_tides_star2 = data->k_div_T_tides_star2;    
    double k_div_T_tides_star3 = data->k_div_T_tides_star3;

    double m1dot = data->m1dot;
    double gamma1 = data->gamma1;
    double gamma2 = data->gamma2;
    double mu1 = data->mu1;
    double mu2 = data->mu2;
    
    double f_quad = data->f_quad;
    double f_oct = data->f_oct;
    double f_tides = data->f_tides;
    double f_mass_transfer = data->f_mass_transfer;
    double f_1PN_1 = data->f_1PN_1;
    double f_1PN_2 = data->f_1PN_2;
    double f_1PN_12 = data->f_1PN_12;
    double f_25PN_1 = data->f_25PN_1;
    double f_25PN_2 = data->f_25PN_2;

	/*	the ODE variables	*/
	double x = Ith(yev,1);
	double y = Ith(yev,2);
	double g1 = Ith(yev,3);
	double g2 = Ith(yev,4);
	double a1 = Ith(yev,5);
	double a2 = Ith(yev,6);
	double cositot = Ith(yev,7);
	double stellarspin1 = Ith(yev,8);
	double stellarspin2 = Ith(yev,9);

	double e1 = 1.0 - pow(10.0,x);
	double e2 = 1.0 - pow(10.0,y);

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
    
	/*	eccentricity functions	*/
	double e1_p2 = e1*e1;
	double e1_p4 = e1_p2*e1_p2;
	double l1_p2 = 1.0 - e1_p2;
    double l1 = sqrt(l1_p2);
    double l1_p3 = l1_p2*l1;

	double e2_p2 = e2*e2;
	double e2_p4 = e2_p2*e2_p2;
	double l2_p2 = 1.0 - e2_p2;
    double l2 = sqrt(l2_p2);
    double l2_p3 = l2_p2*l2;

	double fme1 = (2.0 - 5.0*e1_p2)*(m1*m1 + m2*m2) - 3.0*(2.0 - e1_p2)*m1*m2;
	double f_GR_adot1 = 1.0 + (73.0/24.0)*e1_p2 + (37.0/96.0)*e1_p4;
	double f_GR_edot1 = 1.0 + (121.0/304.0)*e1_p2;

	double f_GR_adot2 = 1.0 + (73.0/24.0)*e2_p2 + (37.0/96.0)*e2_p4;
	double f_GR_edot2 = 1.0 + (121.0/304.0)*e2_p2;
	
    /* tides quantities */
	double R1_div_a1 = R1/a1;
	double R1_div_a1_p2 = R1_div_a1*R1_div_a1;
	double R1_div_a1_p5 = R1_div_a1_p2*R1_div_a1_p2*R1_div_a1;
	double R1_div_a1_p6 = R1_div_a1*R1_div_a1_p5;
	double R2_div_a1 = R2/a1;
	double R2_div_a1_p2 = R2_div_a1*R2_div_a1;
	double R2_div_a1_p5 = R2_div_a1_p2*R2_div_a1_p2*R2_div_a1;
	double R2_div_a1_p6 = R2_div_a1*R2_div_a1_p5;
	
	double tcqr1 = m2_div_m1*R1_div_a1_p6*k_div_T_tides_star1;		/* Common combination of quantities in Hut's equations */
	double tcqr2 = m1_div_m2*R2_div_a1_p6*k_div_T_tides_star2;
    
    double a1_p2 = a1*a1;
    double a1_p3 = a1_p2*a1;
	double n1 = CONST_G*sqrt(m1_plus_m2/a1_p3); // mean orbital angular speed
    double stellarspin1_div_n1 = stellarspin1/n1;
    double stellarspin2_div_n1 = stellarspin2/n1;

    double f_tides1_1 = 
//	double f_tides1 = 1.0 + e1_p2*(c_31div2 + e1_p2*(c_255div8 + e1_p2*(c_185_div16 + e1_p2*(
    
//     + c_255div8*e1_f4 + (185.0/16.0)*e1_f6 + (25.0/64.0)*e1_f8;
//	double f_tides2 = 1.0 + (15.0/2.0)*e1_f2 + (45.0/8.0)*e1_f4 + (5.0/16.0)*e1_f6;
//	double f_tides3 = 1.0 + (15.0/4.0)*e1_f2 + (15.0/8.0)*e1_f4 + (5.0/64.0)*e1_f6;
//	double f_tides4 = 1.0 + (3.0/2.0)*e1_f2 + (1.0/8.0)*e1_f4;
//	double f_tides5 = 1.0 + 3.0*e1_f2 + (3.0/8.0)*e1_f4;

    /* mass transfer quantities */
    double r_A1_P_div_a1,r_A2_P_div_a1,cos_Phi_P,X_L1_P_approx;
    if (f_mass_transfer != 0.0)
    {
        double f=0.0;
        X_L1_P_approx = 0.529 + 0.231*log(m1_div_m2) - f*f*(0.031 + 0.025*e1)*(1.0 + 0.4*log(m1_div_m2)); // 2007ApJ...667.1170S Eq. A15
        r_A1_P_div_a1 = X_L1_P_approx*(1.0-e1);
        
        r_A2_P_div_a1 = R2/a1;
        r_A2_P_div_a1 = 0.0;
        cos_Phi_P = 1.0;        
    }
    
	/*	triple quantities */
	double L1 = m1*m2*sqrt(CONST_G*a1/(m1+m2));
	double L2 = (m1+m2)*m3*sqrt(CONST_G*a2/(m1+m2+m3));
	double G1 = L1*sqrt(l1_p2);
	double G2 = L2*sqrt(l2_p2);

	double a1_div_a2 = a1/a2;
	double C2 = CONST_G*(1.0/16.0)*(m_prod_triple/m_tot_inner)*pow(l2_p2,-3.0/2.0)*a1_div_a2*a1_div_a2/a2;
	double C3 = CONST_G*(-15.0/16.0)*(1.0/4.0)*(m_prod_triple/(m_tot_inner*m_tot_inner))*(m1-m2)*pow(l2_p2,-5.0/2.0)*a1_div_a2*a1_div_a2*a1_div_a2/a2;
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

	double sing1 = sin(g1);
	double sing1d = sin(2.0*g1);	/* 'd' stands for double angle */
	double sing2 = sin(g2);
	double cosg1 = cos(g1);
	double cosg1d = cos(2.0*g1);
	double cosg2 = cos(g2);
	
	/*	Required for octupole order terms	*/
	double B = 2.0 + 5.0*e1_p2 - 7.0*e1_p2*cosg1d;
	double A = 4.0 + 3.0*e1_p2 - (5.0/2.0)*B*sinitot_p2;
	double cosphi = -cosg1*cosg2 - cositot*sing1*sing2;
	

	/* ================================================================
	* The calculations of the ODE right hand sides
	* References: Ford et al. (2000); Blaes et al. (2002); Naoz et al. (2012)
	*/
	

    /*******************************
     * e1dot                       *
     * *****************************/
    
    /* Newtonian point particle -- up and including octupole order */
	double e1dot_newtonian = C2*(l1_p2/G1)*(30.0*e1*sinitot_p2*sing1d) \
		+ C3*e2*(l1_p2/G1)*(35.0*cosphi*sinitot_p2*e1_p2*sing1d \
			- 10.0*cositot*sinitot_p2*cosg1*sing2*l1_p2 \
			- A*(sing1*cosg2 - cositot*cosg1*sing2));
    
    /* post-Newtonian point particle */
    double e1dot_GR_1PN_12=0.0,e1dot_GR_25PN_1=0.0;
    if (f_1PN_12 != 0.0)
    {
        e1dot_GR_1PN_12 = (-9.0/16.0)*(CONST_G*CONST_G*pow(CONST_C_LIGHT,-2.0))*a1*e1*sqrt(l1_p2)*m1*m2*(m1*m1 + m1*m2 + m2*m2) \
            *m3*sinitot_p2*sing1d*pow( pow(a2,3.0)*pow(l2_p2,3.0/2.0)*L1*pow(m1+m2,2.0),-1.0);
    }
    if (f_25PN_1 != 0.0)
    {
        e1dot_GR_25PN_1 = (-304.0/15.0)*pow(CONST_G,3.0)*m_prod_tot_inner*e1*pow(a1,-4.0)*pow(CONST_C_LIGHT,-5.0)*pow(l1_p2,-5.0/2.0)*f_GR_edot1;
    }
    
    /* tides */
    double e1dot_tides_star1=0,e1dot_tides_star2=0;
    if (f_tides != 0.0)
    {
//        e1dot_tides_star1 = -27.0*(1.0+m2_div_m1)*tcqr1*R1_div_a1_p2*e1*pow(l1,-13.0)*(f3 \
//            - (11.0/18.0)*l1_p3*f4*stellarspin1_div_n1);
//        e1dot_tides_star2 = -27.0*(1.0+m1_div_m2)*tcqr2_f*ra2_f2*e1_f*pow(e1_f2com,-13.0/2.0)*(f3_f \
//            - (11.0/18.0)*pow(e1_f2com,3.0/2.0)*f4_f*omegaquot2_f);
    }
    
    /* mass transfer */
    double e1dot_mass_transfer=0.0;
    if (f_mass_transfer != 0.0)
    {
        e1dot_mass_transfer = (l1/(2.0*M_PI))*(m1dot/m1)*( gamma1*m1_div_m2*r_A2_P_div_a1*cos_Phi_P + r_A1_P_div_a1 \
            + 2.0*(gamma1*m1_div_m2-1.0)*(1.0-e1) + 2.0*(1.0-gamma1)*(mu1 + c_1div2)*(1.0-e1)*m1_div_m2/(1.0+m1_div_m2) );
//        if (e1 <= 0.0)
//        {
//            e1dot_mass_transfer = 0.0;
//        }
    }

    /* combined */
	double e1dot = e1dot_newtonian + f_1PN_12*e1dot_GR_1PN_12 + f_25PN_1*e1dot_GR_25PN_1 + f_mass_transfer*e1dot_mass_transfer;
	Ith(ydot,1) = -1.0*pow(10.0,-x)*e1dot/log(10.0);


    /*******************************
     * e2dot                       *
     * *****************************/

    /* Newtonian point particle -- up and including octupole order */
	double e2dot_newtonian = -C3*e1*(l2_p2/G2)*(10.0*cositot*sinitot_p2*l1_p2*sing1*cosg2 \
		+ A*(cosg1*sing2 - cositot*sing1*cosg2));
        
    /* post-Newtonian point particle */        
	double e2dot_GR_25PN_2=0.0;
    if (f_25PN_2 != 0.0)
    {
        f_25PN_2 = (-304.0/15.0)*pow(CONST_G,3.0)*m_tot_inner*m3*m_tot_triple*e2*pow(a2,-4.0)*pow(CONST_C_LIGHT,-5.0)*pow(l2_p2,-5.0/2.0)*f_GR_edot2;
    }
    
    /* combined */
	double e2dot = e2dot_newtonian + f_25PN_2*e2dot_GR_25PN_2;
	Ith(ydot,2) = -1.0*pow(10.0,-y)*e2dot/log(10.0);


    /*******************************
     * g1dot                       *
     * *****************************/
     
    /* Newtonian point particle -- up and including octupole order */
	double g1dot_newtonian = 6.0*C2*((1.0/G1)*(4.0*cositot_p2 + (5.0*cosg1d - 1.0)*(l1_p2 - cositot_p2)) \
			+ (cositot/G2)*(2.0 + e1_p2*(3.0 - 5.0*cosg1d))) \
		- C3*e2*(e1*((1.0/G2) + (cositot/G1))*(sing1*sing2*(10.0*(3.0*cositot_p2 \
			- 1.0)*(1.0 - e1_p2) + A) - 5.0*B*cositot*cosphi) \
			- (l1_p2/(e1*G1))*(sing1*sing2*10.0*cositot*sinitot_p2*(1.0 - 3.0*e1_p2) \
			+ cosphi*(3.0*A - 10.0*cositot_p2 + 2.0)));

    /* post-Newtonian point particle */  
    double g1dot_GR_1PN_1=0.0,g1dot_GR_1PN_12=0.0;
    if (f_1PN_1 != 0.0)
    {
        g1dot_GR_1PN_1 = 3.0*pow(CONST_C_LIGHT_P2*a1*l1_p2,-1.0)*pow(CONST_G*m_tot_inner/a1,3.0/2.0);
    }
    if (f_1PN_12 != 0.0)
    {
        g1dot_GR_1PN_12 = (CONST_G*pow(CONST_C_LIGHT,-2.0))*pow( 4.0*pow(a2,3.0)*pow(l2_p2,3.0/2.0)*(m1+m2),-1.0) \
            *( G2*( 8.0*(m1+m2) + 6.0*m3)*cositot \
                - (l1_p2/G1)*(CONST_G)*a1*m1*m2*m3*pow(8.0*(m1+m2),-1.0)*(tildefme1*(1.0 - 3.0*cositot_p2) \
                    + 18.0*(m1*m1 + m1*m2 + m2*m2)*cosg1d*sinitot_p2) \
                - ( (cositot/G1) + (1.0/G2) )*( G1*G2*(8.0*(m1+m2) + 6.0*m3) \
                    + cositot*a1*CONST_G*m1*m2*m3*pow(8.0*(m1+m2),-1.0) \
                    *(6.0*fme1 + 18.0*e1_p2*(m1*m1 + m1*m2 + m2*m2)*cosg1d) ) );
    }
	
    /* mass transfer */    
    double g1dot_mass_transfer=0.0;
    if (f_mass_transfer != 0.0)
    {
        g1dot_mass_transfer = 0.0;
    }
    
    /* combined */    
	Ith(ydot,3) = g1dot_newtonian + f_1PN_1*g1dot_GR_1PN_1 + f_1PN_12*g1dot_GR_1PN_12 + f_mass_transfer*g1dot_mass_transfer;


    /*******************************
     * g2dot                       *
     * *****************************/
     
    /* Newtonian point particle -- up and including octupole order */
	double g2dot_newtonian = 3.0*C2*((2.0*cositot/G1)*(2.0 + e1_p2*(3.0 - 5.0*cosg1d)) \
			+ (1.0/G2)*(4.0 + 6.0*e1_p2 + (5.0*cositot_p2 - 3.0)*(2.0 + e1_p2*(3.0 - 5.0*cosg1d)))) \
		+ C3*e1*(sing1*sing2*(((4.0*e2_p2 + 1.0)/(e2*G2))*10.0*cositot*sinitot_p2*l1_p2 \
			- e2*((1.0/G1) + (cositot/G2))*(A + 10.0*(3.0*cositot_p2 - 1.0)*l1_p2)) \
			+ cosphi*(5.0*B*cositot*e2*((1.0/G1) + (cositot/G2)) + ((4.0*e2_p2 + 1.0)/(e2*G2))*A));

    /* post-Newtonian point particle */  
    double g2dot_GR_1PN_2=0.0,g2dot_GR_1PN_12=0.0;
    if (f_1PN_1 != 0.0)
    {    
        g2dot_GR_1PN_2 = 3.0*pow(CONST_C_LIGHT_P2*a2*l2_p2,-1.0)*pow(CONST_G*m_tot_triple/a2,3.0/2.0);
    }
    if (f_1PN_12 != 0.0)
    {
        g2dot_GR_1PN_12 = (CONST_G*pow(CONST_C_LIGHT_P2,-1.0))*pow( 4.0*pow(a2,3.0)*pow(l2_p2,3.0/2.0)*(m1+m2),-1.0) \
            *( -2.0*G1*(8.0*(m1+m2)+6.0*m3)*cositot - 3.0*CONST_G*a1*m1*m2*m3*pow(8.0*(m1+m2)*G2,-1.0) \
                    *( fme1 - 3.0*fme1*cositot_p2 + 9.0*e1_p2*(m1*m1+m1*m2+m2*m2)*cosg1d*sinitot_p2 ) \
                + ( (cositot/G2) + (1.0/G1) )*(G1*G2*(8.0*(m1+m2)+6.0*m3) \
                    + CONST_G*a1*m1*m2*m3*pow(8.0*(m1+m2),-1.0)*cositot*(6.0*fme1 \
                        + 18.0*e1_p2*(m1*m1+m1*m2+m2*m2)*cosg1d)) );
    }
    
    /* combined */
	Ith(ydot,4) = g2dot_newtonian + f_1PN_2*g2dot_GR_1PN_2 + f_1PN_12*g2dot_GR_1PN_12;		
	
   
    /*******************************
     * a1dot                       *
     * *****************************/

    /* post-Newtonian point particle */  
	double a1dot_GR_25PN_1=0.0;
    if (f_25PN_1 != 0.0)
    {
        a1dot_GR_25PN_1 = (-64.0/5.0)*pow(CONST_G,3.0)*pow(CONST_C_LIGHT,-5.0)*m_prod_tot_inner*pow(a1,-3.0)*pow(l1_p2,-7.0/2.0)*f_GR_adot1;
    }

    /* mass transfer */
    double a1dot_mass_transfer=0.0;
    if (f_mass_transfer != 0.0)
    {
        a1dot_mass_transfer = (a1/M_PI)*(m1dot/m1)*(1.0/l1)*( e1*r_A1_P_div_a1 + gamma1*m1_div_m2*e1*r_A2_P_div_a1*cos_Phi_P \
            + (gamma1*m1_div_m2-1.0)*(1.0 - e1_p2) \
            + (1.0-gamma1)*(mu1 + c_1div2)*l1_p2*m1_div_m2/(1.0 + m1_div_m2) );
    }
        
    /* combined */
	Ith(ydot,5) = f_25PN_1*a1dot_GR_25PN_1 + f_mass_transfer*a1dot_mass_transfer;


    /*******************************
     * a2dot                       *
     * *****************************/

    /* post-Newtonian point particle */  
	double a2dot_GR_25PN_2=0.0;
    if (f_25PN_2 != 0.0)
    {
        a2dot_GR_25PN_2 = (-64.0/5.0)*pow(CONST_G,3.0)*pow(CONST_C_LIGHT,-5.0)*m_tot_inner*m3*m_tot_triple*pow(a2,-3.0)*pow(l2_p2,-7.0/2.0)*f_GR_adot2;
    }

    /* combined */
	Ith(ydot,6) = f_25PN_2*a2dot_GR_25PN_2;


    /************************************
     * cositotdot                       *
     * due to triple interaction only!  *
     * **********************************/

	double G1dot = (-1.0)*G1*e1*(e1dot_newtonian+f_1PN_12*e1dot_GR_1PN_12)/l1_p2;
	double G2dot = (-1.0)*G2*e2*e2dot_newtonian/l2_p2;
	double cositotdot = (-1.0/(G1*G2))*(G1dot*(G1 + G2*cositot) + G2dot*(G2 + G1*cositot));

	Ith(ydot,7) = cositotdot;


    /************************************
     * stellarspin1dot                  *
     * **********************************/

	Ith(ydot,8) = 0.0;


    /************************************
     * stellarspin2dot                  *
     * **********************************/

	Ith(ydot,9) = 0.0;

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
