/*	Worker code for SecularTriple, a secular triple gravitational dynamics code taking into account Newtonian, 1PN and 2.5PN terms	*/
/*	The relevant ODEs are solved consistently for each user supplied timestep using CVODE (Cohen & Hindmarsh 1996)	*/

#include <math.h>
#include <cstdlib>
#include "code.h"

//#include "cvode/cvode.h"				    /* prototypes for CVODE fcts., consts. */
//#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
//#include "cvode/cvode_dense.h"			    /* prototype for CVDense */
//#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
//#include "cvode/sundials_types.h"			/* definition of type realtype */

#define CONST_G			    (double)	6.67259e-8
#define CONST_C_LIGHT		(double)	2.99792458e10
#define CONST_C_LIGHT_POW2	(double)	CONST_C_LIGHT*CONST_C_LIGHT

/*	ODE solver related quantities	*/
#define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

//#define NEQ		    9				        /* number of ODE equations */
#define RTOL		RCONST(1.0e-10)			/* scalar relative tolerance - acceptable value is 1.0e-10 as determined by trial and error */
#define ATOL1		RCONST(1.0e-8)
#define ATOL2		RCONST(1.0e-8)
#define ATOL3		RCONST(1.0e-8)
#define ATOL4		RCONST(1.0e-8)
#define ATOL5		RCONST(1.0e1)
#define ATOL6		RCONST(1.0e-8)
#define ATOL7		RCONST(1.0e-8)
#define ATOL8		RCONST(1.0e-8)	
#define ATOL9		RCONST(1.0e-8)				
#define TIMESTEP0	RCONST(1.0)				/* initial internal ODE timestep */
#define MAXNUMSTEPS 	5e8					/* maximum number of internal steps */
#define MAXTIME		RCONST(13.7e9*365.25*24.0*3600.0)	/* maximum integration time */
#define MAXNUMCONVFAIL	20					/* maximum number of convergence failures */

typedef struct {
	double m1,m2,m3,R1,R2,R3;
} *UserData;

static int fev(realtype t, N_Vector yev, N_Vector ydot, void *data);
static int froot(realtype t, N_Vector yev, realtype *gout, void *data);	/* ODE Rootfinding function */
double f_oct = 1.0;
int equations_of_motion_specification = 0;
double f_1PN_1 = 0;
double f_1PN_2 = 0;
double f_1PN_12 = 0;
double f_25PN_1 = 0;
double f_25PN_2 = 0;
double relative_tolerance = 1.0e-10;

/*	interface parameters	*/
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

int get_f_oct(double *f_oct_t)
{
	*f_oct_t = f_oct;
	return 0;
}
int set_f_oct(double f_oct_t)
{
	f_oct = f_oct_t;
	return 0;
}

int get_f_1PN_1(double *f_1PN_1_t)
{
	*f_1PN_1_t = f_1PN_1;
	return 0;
}
int set_f_1PN_1(double f_1PN_1_t)
{
	f_1PN_1 = f_1PN_1_t;
	return 0;
}

int get_f_1PN_2(double *f_1PN_2_t)
{
	*f_1PN_2_t = f_1PN_2;
	return 0;
}
int set_f_1PN_2(double f_1PN_2_t)
{
	f_1PN_2 = f_1PN_2_t;
	return 0;
}

int get_f_1PN_12(double *f_1PN_12_t)
{
	*f_1PN_12_t = f_1PN_12;
	return 0;
}
int set_f_1PN_12(double f_1PN_12_t)
{
	f_1PN_12 = f_1PN_12_t;
	return 0;
}

int get_f_25PN_1(double *f_25PN_1_t)
{
	*f_25PN_1_t = f_25PN_1;
	return 0;
}
int set_f_25PN_1(double f_25PN_1_t)
{
	f_25PN_1 = f_25PN_1_t;
	return 0;
}

int get_f_25PN_2(double *f_25PN_2_t)
{
	*f_25PN_2_t = f_25PN_2;
	return 0;
}
int set_f_25PN_2(double f_25PN_2_t)
{
	f_25PN_2 = f_25PN_2_t;
	return 0;
}


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
    int * out_flag, int * error_flag)
{

    if (INCL1==0.0) { INCL1 = 1.0e-15; }
    if (INCL2==0.0) { INCL2 = 1.0e-15; }
    if (AP1==0.0) { AP1 = 1.0e-15; }
    if (AP2==0.0) { AP2 = 1.0e-15; }
    if (LAN1==0.0) { LAN1 = 1.0e-15; }
    if (LAN2==0.0) { LAN2 = 1.0e-15; }

    /*********************************************************************
     * quantities that are assumed to be constant during the integration *
     *********************************************************************/
	UserData data;
	data = NULL;
	data = (UserData) malloc(sizeof *data);
	data->m1 = m1;
	data->m2 = m2;
	data->m3 = m3;
	data->R1 = R1;
	data->R2 = R2;
    data->R3 = R3;

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
    double e1_vec[3], e2_vec[3], h1_vec[3], h2_vec[3], q1_vec[3], q2_vec[3];    
    double cos_INCL1 = cos(INCL1);
    double cos_INCL2 = cos(INCL2);
    double sin_INCL1 = sin(INCL1);
    double sin_INCL2 = sin(INCL2);
    double cos_INCL12 = sin_INCL1*sin_INCL2*cos(LAN1-LAN2) + cos_INCL1*cos_INCL2;
    double INCL12 = acos(cos_INCL12);

    if (equations_of_motion_specification==0)
    {
        /********************************************************
         * solve equations of motion based on delaunay elements *
         * ******************************************************/

        NEQ = 9;
        yev = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev, "N_VNew_Serial", 0)) return;
        yev_out = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev_reached, "N_VNew_Serial", 0)) return;
        abstol = N_VNew_Serial(NEQ); 
//	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return;

        Ith(yev,1) = log10(1.0 - e1);
        Ith(yev,2) = log10(1.0 - e2);
        Ith(yev,3) = AP1;		/* Inner orbit argument of periastron - unit: rad */
        Ith(yev,4) = AP2;		/* Outer orbit argument of periastron - unit: rad */
        Ith(yev,5) = a1;		/* Inner orbit semi-major axis - unit: cm */
        Ith(yev,6) = a2;		/* Outer orbit semi-major axis - unit: cm */
        Ith(yev,7) = cos_INCL12;
        Ith(yev,8) = 0.0;	/* Spin angular frequency of star 1 - unit: rad s^-1 */
        Ith(yev,9) = 0.0;	/* Spin angular frequency of star 2 - unit: rad s^-1 */

        Ith(abstol,1) = ATOL1;   
        Ith(abstol,2) = ATOL2;
        Ith(abstol,3) = ATOL3;
        Ith(abstol,4) = ATOL4;
        Ith(abstol,5) = ATOL5;   
        Ith(abstol,6) = ATOL6;
        Ith(abstol,7) = ATOL7;
        Ith(abstol,8) = ATOL8;
        Ith(abstol,9) = ATOL9;
        

        cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        flag = CVodeInit(cvode_mem, fev_delaunay, t, yev);        
        
    }
    else
    {
        /********************************************************************
         * solve equations of motion based on triad of unit vectors (q,h,e) *
         * ******************************************************************/        
         
        NEQ = 12;
        yev = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev, "N_VNew_Serial", 0)) return;
        yev_out = N_VNew_Serial(NEQ);
//	if (check_flag((void *)yev_reached, "N_VNew_Serial", 0)) return;
        abstol = N_VNew_Serial(NEQ); 
//	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return;         

        compute_triad_vectors_from_orbital_elements(m1,m2,m3,a1,a2,e1,e2,INCL1,INCL2,AP1,AP2,LAN1,LAN2,
            e1_vec,e2_vec,h1_vec,h2_vec,q1_vec,q2_vec);
        FILE *fp;
        fp = fopen("/data2/amuse-svn/output_stream.txt","w");
        fprintf(fp,"test %g %g %g ",e1_vec[2],h1_vec[2],q1_vec[2]);
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

	int nroot = 2;
	int rootsfound[nroot];
	flag = CVodeRootInit(cvode_mem, nroot, froot);
//	if (check_flag(&flag, "CVodeRootInit", 1)) return;	

	double t_end = t + dt;
	double t_end_cvode;
	int out_flag_f = 0;
	int error_flag_f = 0;

	flag_s = CVode(cvode_mem, t_end, yev_out, &t_end_cvode, CV_NORMAL);	
	if (flag_s == CV_SUCCESS)
	{
		out_flag_f = 0;
		error_flag_f = 0;
	}
	else if (flag_s == CV_ROOT_RETURN)
	{
		CVodeGetRootInfo(cvode_mem,rootsfound);
		if (rootsfound[0] == 1 || rootsfound[0] == -1)
		{
			out_flag_f = 1;
		}
		if (rootsfound[1] == 1 || rootsfound[1] == -1)
		{
			out_flag_f = 2;
		}
		error_flag_f = 0;
	}
	else
		out_flag_f = 99;
		error_flag_f = flag_s;


    /**********
     * output *
     * ********/
     
	*m1_out = m1;
	*m2_out = m2;
	*m3_out = m3;
	*R1_out = R1;
	*R2_out = R2;
    *R3_out = R3;
    
	*t_out = t_end_cvode;
	*out_flag = out_flag_f;
	*error_flag = error_flag_f;

    double x_out,y_out;
//      double L1_out,L2_out,G1_out,G2_out,G12_out;
    double cos_INCL1_out;
    double cos_INCL2_out;
    double sin_INCL1_out;
    double sin_INCL2_out;
    double cos_INCL12_out;
//    double INCL12 = acos(cos_INCL12);


    if (equations_of_motion_specification==0)
    {
        x_out = Ith(yev_out,1);
        y_out = Ith(yev_out,2);
        
        *e1_out = 1.0 - pow(10,x_out);
        *e2_out = 1.0 - pow(10,y_out);
        *AP1_out = Ith(yev_out,3);
        *AP2_out = Ith(yev_out,4);
        *a1_out = Ith(yev_out,5);
        *a2_out = Ith(yev_out,6);        

        cos_INCL12_out = Ith(yev_out,7);

        if (cos_INCL12_out > 1.0)
        {
            cos_INCL12_out = 2.0 - cos_INCL12_out;
        }
    
        if (cos_INCL12_out < -1.0)
        {
            cos_INCL12_out = -2.0 - cos_INCL12_out;
        }

        *INCL1_out = acos(cos_INCL12_out);
        *INCL2_out = 0.0;
        *INCL12_out = acos(cos_INCL12_out);
        
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
    }

	return 0;
}

int compute_triad_vectors_from_orbital_elements(double m1, double m2, double m3, double a1, double a2, double e1, double e2, double INCL1, double INCL2, double AP1, double AP2, double LAN1, double LAN2,
    double e1_vec[3], double e2_vec[3], double h1_vec[3], double h2_vec[3], double q1_vec[3], double q2_vec[3])
{
    double L1,L2,G1,G2,G12;
    double cos_INCL1 = cos(INCL1);
    double cos_INCL2 = cos(INCL2);
    double sin_INCL1 = sin(INCL1);
    double sin_INCL2 = sin(INCL2);
    double cos_AP1 = cos(AP1);
    double cos_AP2 = cos(AP2);
    double sin_AP1 = sin(AP1);
    double sin_AP2 = sin(AP2);
    double cos_LAN1 = cos(LAN1);
    double cos_LAN2 = cos(LAN2);
    double sin_LAN1 = sin(LAN1);
    double sin_LAN2 = sin(LAN2);
    
    L1 = m1*m2*sqrt(CONST_G*a1/(m1+m2));
    L2 = (m1+m2)*m3*sqrt(CONST_G*a2/(m1+m2+m3));
    G1 = L1*sqrt(1.0 - e1*e1);
    G2 = L2*sqrt(1.0 - e2*e2);
//    G12 = sqrt(G1*G1 + G2*G2 + 2.0*G1*G2*cos(INCL12));
//    cos_INCL1 = (G12*G12 + G1*G1 - G2*G2)/(2.0*G12*G1);
//    cos_INCL2 = (G12*G12 - G1*G1 + G2*G2)/(2.0*G12*G2);
//    INCL1 = acos(cos_INCL1);
//    INCL2 = acos(cos_INCL2);
//    sin_INCL1 = sin(INCL1);
//    sin_INCL2 = sin(INCL2);
//    sin_INCL1 = sqrt(1.0 - cos_INCL1*cos_INCL1);
//    sin_INCL2 = sqrt(1.0 - cos_INCL2*cos_INCL2);
    
    e1_vec[0] = e1*(cos_LAN1*cos_AP1 - sin_LAN1*sin_AP1*cos_INCL1);
    e1_vec[1] = e1*(sin_LAN1*cos_AP1 + cos_LAN1*sin_AP1*cos_INCL1);
    e1_vec[2] = e1*(sin_AP1*sin_INCL1);
    
    h1_vec[0] = G1*(sin_LAN1*sin_INCL1);
    h1_vec[1] = G1*(-cos_LAN1*sin_INCL1);
    h1_vec[2] = G1*(cos_INCL1);

    e2_vec[0] = e2*(cos_LAN2*cos_AP2 - sin_LAN2*sin_AP2*cos_INCL2);
    e2_vec[1] = e2*(sin_LAN2*cos_AP2 + cos_LAN2*sin_AP2*cos_INCL2);
    e2_vec[2] = e2*(sin_AP2*sin_INCL2);
    
    h2_vec[0] = G2*(sin_LAN2*sin_INCL2);
    h2_vec[1] = G2*(-cos_LAN2*sin_INCL2);
    h2_vec[2] = G2*(cos_INCL2);

    double h1_sq = norm3_squared(h1_vec);
    double h2_sq = norm3_squared(h2_vec);
    
    double h1 = sqrt(h1_sq);
    double h2 = sqrt(h2_sq);
    
    double cos_INCL12_temp = dot3(h1_vec,h2_vec)/(h1*h2);
    
    cross3(h1_vec,e1_vec,q1_vec);
    cross3(h2_vec,e2_vec,q2_vec);

    FILE *fp;
    fp = fopen("/data2/amuse-svn/output_stream2.txt","w");
    fprintf(fp,"test %g %g %g %g %g %g",
        pow(e1_vec[0]/e1,2.0) + pow(e1_vec[1]/e1,2.0) + pow(e1_vec[2]/e1,2.0),
        pow(h1_vec[0]/G1,2.0) + pow(h1_vec[1]/G1,2.0) + pow(h1_vec[2]/G1,2.0),
        pow(e2_vec[0]/e2,2.0) + pow(e2_vec[1]/e2,2.0) + pow(e2_vec[2]/e2,2.0),
        pow(h2_vec[0]/G2,2.0) + pow(h2_vec[1]/G2,2.0) + pow(h2_vec[2]/G2,2.0),
        h2_vec[2]);
//        cos(INCL12), sin_INCL1*sin_INCL2*cos(LAN1-LAN2) + cos_INCL1*cos_INCL2);
    fclose(fp);    
    
    return 0;
}

int compute_orbital_elements_from_triad_vectors(double m1, double m2, double m3,
    double e1_vec[3],double e2_vec[3],double h1_vec[3],double h2_vec[3],double q1_vec[3],double q2_vec[3],
    double *a1_out, double *a2_out, double *e1_out, double *e2_out, double *INCL1_out, double *INCL2_out, double *AP1_out, double *AP2_out, double *LAN1_out, double *LAN2_out)
{
    double e1_sq = norm3_squared(e1_vec);
    double e2_sq = norm3_squared(e2_vec);    
    double e1 = sqrt(e1_sq);
    double e2 = sqrt(e2_sq);
    
    double h1_sq = norm3_squared(h1_vec);
    double h2_sq = norm3_squared(h2_vec);
    
    double a1 = h1_sq*(m1+m2)/( CONST_G*m1*m1*m2*m2*(1.0 - e1_sq) );
    double a2 = h2_sq*(m1+m2+m3)/( CONST_G*(m1+m2)*(m1+m2)*m3*m3*(1.0 - e2_sq) );
    
    double h1 = sqrt(h1_sq);
    double h2 = sqrt(h2_sq);
    
    double x_vec[3] = {1.0,0.0,0.0};
    double y_vec[3] = {0.0,1.0,0.0};
    double z_vec[3] = {0.0,0.0,1.0};

    double cos_INCL1 = dot3(h1_vec,z_vec)/h1;
    double cos_INCL2 = dot3(h2_vec,z_vec)/h2;
//    double cos_INCL12 = dot3(h1_vec,h2_vec)/(h1*h2);

    double LAN1_vec[3],LAN2_vec[3];
    double LAN1_vec_unit[3],LAN2_vec_unit[3];
    cross3(z_vec,h1_vec,LAN1_vec);
    cross3(z_vec,h2_vec,LAN2_vec);
    double LAN1_vec_norm = norm3(LAN1_vec);
    double LAN2_vec_norm = norm3(LAN2_vec);

    double e1_vec_unit[3], e2_vec_unit[3];
    double h1_vec_unit[3], h2_vec_unit[3];
    int i=0;
    for (int i=0; i<3; i++)
    {
        LAN1_vec_unit[i] = LAN1_vec[i]/LAN1_vec_norm;
        LAN2_vec_unit[i] = LAN2_vec[i]/LAN2_vec_norm;
        e1_vec_unit[i] = e1_vec[i]/e1;
        e2_vec_unit[i] = e2_vec[i]/e2;
        h1_vec_unit[i] = h1_vec[i]/h1;
        h2_vec_unit[i] = h2_vec[i]/h2;
    }

    double sin_LAN1 = dot3(LAN1_vec,y_vec);
    double sin_LAN2 = dot3(LAN2_vec,y_vec);    
    double cos_LAN1 = dot3(LAN1_vec,x_vec);
    double cos_LAN2 = dot3(LAN2_vec,x_vec);    

    double e1_cross_h1[3], e2_cross_h2[3];
    cross3(e1_vec_unit,h1_vec_unit,e1_cross_h1);
    cross3(e2_vec_unit,h2_vec_unit,e2_cross_h2);
    double sin_AP1 = dot3(LAN1_vec_unit,e1_cross_h1);
    double sin_AP2 = dot3(LAN2_vec_unit,e2_cross_h2);
    double cos_AP1 = dot3(LAN1_vec_unit,e1_vec_unit);
    double cos_AP2 = dot3(LAN2_vec_unit,e2_vec_unit);

    double LAN1 = atan2(sin_LAN1,cos_LAN1);
    double LAN2 = atan2(sin_LAN2,cos_LAN2);    
    double AP1 = atan2(sin_AP1,cos_AP1);
    double AP2 = atan2(sin_AP2,cos_AP2);

    *a1_out = a1;
    *a2_out = a2;
    *e1_out = e1;
    *e2_out = e2;
    *INCL1_out = acos(cos_INCL1);
    *INCL2_out = acos(cos_INCL2);
    *AP1_out = AP1;
    *AP2_out = AP2;
    *LAN1_out = LAN1;
    *LAN2_out = LAN2;
    
    return 0;
}

void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}


static int fev_triad(realtype t, N_Vector yev, N_Vector ydot, void *data)
{
	
	UserData data_f;
	data_f = (UserData) data;
	
	/*	Constants which appear in the ODE right hand sides	*/
	double m1_f = data_f->m1;					
	double m2_f = data_f->m2;				
	double m3_f = data_f->m3;				

	/*	The ODE variables	*/
	double x_f = Ith(yev,1);
	double y_f = Ith(yev,2);
	double g1_f = Ith(yev,3);
	double g2_f = Ith(yev,4);
	double a1_f = Ith(yev,5);
	double a2_f = Ith(yev,6);
	double cositot_f = Ith(yev,7);
	double omegaspin1_f = Ith(yev,8);
	double omegaspin2_f = Ith(yev,9);

	double e1_f = 1.0 - pow(10.0,x_f);
	double e2_f = 1.0 - pow(10.0,y_f);

	/*	Some mass combinations	*/
	double q_f = m2_f/m1_f;
	double qinv_f = 1.0/q_f;
	double m_prod_inner_f = m1_f*m2_f;
	double m_tot_inner_f = m1_f+m2_f;
	double m_prod_triple_f = m_prod_inner_f*m3_f;
	double m_tot_triple_f = m_tot_inner_f+m3_f;
	double m_prod_tot_inner_f = m_prod_inner_f*m_tot_inner_f;
	double tildefme1_f = -10.0*m1_f*m1_f + 6.0*m1_f*m2_f - 10.0*m2_f*m2_f;

	/*	Eccentricity functions	*/
	double e1_f2 = e1_f*e1_f;
	double e1_f4 = e1_f2*e1_f2;
	double e1_f2com = 1.0 - e1_f2;		/* 'com' stands for complement */

	double e2_f2 = e2_f*e2_f;
	double e2_f4 = e2_f2*e2_f2;
	double e2_f2com = 1.0 - e2_f2;

	double fme1_f = (2.0 - 5.0*e1_f2)*(m1_f*m1_f + m2_f*m2_f) - 3.0*(2.0 - e1_f2)*m1_f*m2_f;
	double f_GR_adot1 = 1.0 + (73.0/24.0)*e1_f2 + (37.0/96.0)*e1_f4;
	double f_GR_edot1 = 1.0 + (121.0/304.0)*e1_f2;

	double f_GR_adot2 = 1.0 + (73.0/24.0)*e2_f2 + (37.0/96.0)*e2_f4;
	double f_GR_edot2 = 1.0 + (121.0/304.0)*e2_f2;
}	

static int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data)
{
	
	UserData data_f;
	data_f = (UserData) data;
	
	/*	Constants which appear in the ODE right hand sides	*/
	double m1_f = data_f->m1;					
	double m2_f = data_f->m2;				
	double m3_f = data_f->m3;				

	/*	The ODE variables	*/
	double x_f = Ith(yev,1);
	double y_f = Ith(yev,2);
	double g1_f = Ith(yev,3);
	double g2_f = Ith(yev,4);
	double a1_f = Ith(yev,5);
	double a2_f = Ith(yev,6);
	double cositot_f = Ith(yev,7);
	double omegaspin1_f = Ith(yev,8);
	double omegaspin2_f = Ith(yev,9);

	double e1_f = 1.0 - pow(10.0,x_f);
	double e2_f = 1.0 - pow(10.0,y_f);

	/*	Some mass combinations	*/
	double q_f = m2_f/m1_f;
	double qinv_f = 1.0/q_f;
	double m_prod_inner_f = m1_f*m2_f;
	double m_tot_inner_f = m1_f+m2_f;
	double m_prod_triple_f = m_prod_inner_f*m3_f;
	double m_tot_triple_f = m_tot_inner_f+m3_f;
	double m_prod_tot_inner_f = m_prod_inner_f*m_tot_inner_f;
	double tildefme1_f = -10.0*m1_f*m1_f + 6.0*m1_f*m2_f - 10.0*m2_f*m2_f;

	/*	Eccentricity functions	*/
	double e1_f2 = e1_f*e1_f;
	double e1_f4 = e1_f2*e1_f2;
	double e1_f2com = 1.0 - e1_f2;		/* 'com' stands for complement */

	double e2_f2 = e2_f*e2_f;
	double e2_f4 = e2_f2*e2_f2;
	double e2_f2com = 1.0 - e2_f2;

	double fme1_f = (2.0 - 5.0*e1_f2)*(m1_f*m1_f + m2_f*m2_f) - 3.0*(2.0 - e1_f2)*m1_f*m2_f;
	double f_GR_adot1 = 1.0 + (73.0/24.0)*e1_f2 + (37.0/96.0)*e1_f4;
	double f_GR_edot1 = 1.0 + (121.0/304.0)*e1_f2;

	double f_GR_adot2 = 1.0 + (73.0/24.0)*e2_f2 + (37.0/96.0)*e2_f4;
	double f_GR_edot2 = 1.0 + (121.0/304.0)*e2_f2;
	
	/*	General quantities */
	double L1_f = m1_f*m2_f*sqrt(CONST_G*a1_f/(m1_f+m2_f));
	double L2_f = (m1_f+m2_f)*m3_f*sqrt(CONST_G*a2_f/(m1_f+m2_f+m3_f));
	double Ga1_f = L1_f*sqrt(1.0 - e1_f2);
	double Ga2_f = L2_f*sqrt(1.0 - e2_f2);

	double a1a2quot_f = a1_f/a2_f;
	double C2_f = CONST_G*(1.0/16.0)*m_prod_triple_f*pow(m_tot_inner_f,-1.0)*pow(e2_f2com,-3.0/2.0)*pow(a1a2quot_f,2.0)*pow(a2_f,-1.0);
	double C3_f = CONST_G*(-15.0/16.0)*(1.0/4.0)*m_prod_triple_f*pow(m_tot_inner_f,-2.0)*(m1_f-m2_f)*pow(e2_f2com,-5.0/2.0)*pow(a1a2quot_f,3.0)*pow(a2_f,-1.0);
  C3_f *= f_oct;

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

	double sing1_f = sin(g1_f);
	double sing1_fd = sin(2.0*g1_f);	/* 'd' stands for double angle */
	double sing2_f = sin(g2_f);
	double cosg1_f = cos(g1_f);
	double cosg1_fd = cos(2.0*g1_f);
	double cosg2_f = cos(g2_f);
	
	/*	Required for octupole order terms	*/
	double B_f = 2.0 + 5.0*e1_f2 - 7.0*e1_f2*cosg1_fd;
	double A_f = 4.0 + 3.0*e1_f2 - (5.0/2.0)*B_f*sinitot_f2;
	double cosphi_f = -cosg1_f*cosg2_f - cositot_f*sing1_f*sing2_f;
	

	/* ================================================================
	* The calculations of the ODE right hand sides
	* References: Ford et al. (2000); Blaes et al. (2002); Naoz et al. (2012)
	*/
	

	/* ----------------------------------------------------------------
	* e1dot
	*/ 

	double e1dot_newtonian = C2_f*((1.0-e1_f2)/Ga1_f)*(30.0*e1_f*sinitot_f2*sing1_fd) \
		+ C3_f*e2_f*((1.0-e1_f2)/Ga1_f)*(35.0*cosphi_f*sinitot_f2*e1_f2*sing1_fd \
			- 10.0*cositot_f*sinitot_f2*cosg1_f*sing2_f*(1.0-e1_f2) \
			- A_f*(sing1_f*cosg2_f - cositot_f*cosg1_f*sing2_f));
	double e1dot_GR_1PN_12 = (-9.0/16.0)*(CONST_G*CONST_G*pow(CONST_C_LIGHT,-2.0))*a1_f*e1_f*sqrt(e1_f2com)*m1_f*m2_f*(m1_f*m1_f + m1_f*m2_f + m2_f*m2_f) \
		*m3_f*sinitot_f2*sing1_fd*pow( pow(a2_f,3.0)*pow(e2_f2com,3.0/2.0)*L1_f*pow(m1_f+m2_f,2.0),-1.0);
	double e1dot_GR_25PN_1 = (-304.0/15.0)*pow(CONST_G,3.0)*m_prod_tot_inner_f*e1_f*pow(a1_f,-4.0)*pow(CONST_C_LIGHT,-5.0)*pow(e1_f2com,-5.0/2.0)*f_GR_edot1;

	double e1dot_f = e1dot_newtonian + f_1PN_12*e1dot_GR_1PN_12 + f_25PN_1*e1dot_GR_25PN_1;
	Ith(ydot,1) = -1.0*pow(10.0,-1.0*x_f)*e1dot_f/log(10.0);


	/* ----------------------------------------------------------------
	* e2dot
	*/ 

	double e2dot_newtonian = -C3_f*e1_f*((1.0-e2_f2)/Ga2_f)*(10.0*cositot_f*sinitot_f2*(1.0-e1_f2)*sing1_f*cosg2_f \
		+ A_f*(cosg1_f*sing2_f - cositot_f*sing1_f*cosg2_f));
	double e2dot_GR_25PN_2 = (-304.0/15.0)*pow(CONST_G,3.0)*m_tot_inner_f*m3_f*m_tot_triple_f*e2_f*pow(a2_f,-4.0)*pow(CONST_C_LIGHT,-5.0)*pow(e2_f2com,-5.0/2.0)*f_GR_edot2;
	double e2dot_f = e2dot_newtonian + f_25PN_2*e2dot_GR_25PN_2;
	Ith(ydot,2) = -1.0*pow(10.0,-1.0*y_f)*e2dot_f/log(10.0);


	/* ----------------------------------------------------------------
	* g1dot
	*/ 

	double g1dot_newtonian = 6.0*C2_f*((1.0/Ga1_f)*(4.0*cositot_f2 + (5.0*cosg1_fd - 1.0)*(1.0 - e1_f2 - cositot_f2)) \
			+ (cositot_f/Ga2_f)*(2.0 + e1_f2*(3.0 - 5.0*cosg1_fd))) \
		- C3_f*e2_f*(e1_f*((1.0/Ga2_f) + (cositot_f/Ga1_f))*	(sing1_f*sing2_f*(10.0*(3.0*cositot_f2 \
			- 1.0)*(1.0 - e1_f2) + A_f) - 5.0*B_f*cositot_f*cosphi_f) \
			- ((1.0-e1_f2)/(e1_f*Ga1_f))*(sing1_f*sing2_f*10.0*cositot_f*sinitot_f2*(1.0 - 3.0*e1_f2) \
			+ cosphi_f*(3.0*A_f - 10.0*cositot_f2 + 2.0)));

	double g1dot_GR_1PN_1 = 3.0*pow(CONST_C_LIGHT_POW2*a1_f*e1_f2com,-1.0)*pow(CONST_G*m_tot_inner_f/a1_f,3.0/2.0);
	double g1dot_GR_1PN_12 = (CONST_G*pow(CONST_C_LIGHT,-2.0))*pow( 4.0*pow(a2_f,3.0)*pow(e2_f2com,3.0/2.0)*(m1_f+m2_f),-1.0) \
		*( Ga2_f*( 8.0*(m1_f+m2_f) + 6.0*m3_f)*cositot_f \
			- (e1_f2com/Ga1_f)*(CONST_G)*a1_f*m1_f*m2_f*m3_f*pow(8.0*(m1_f+m2_f),-1.0)*(tildefme1_f*(1.0 - 3.0*cositot_f2) \
				+ 18.0*(m1_f*m1_f + m1_f*m2_f + m2_f*m2_f)*cosg1_fd*sinitot_f2) \
			- ( (cositot_f/Ga1_f) + (1.0/Ga2_f) )*( Ga1_f*Ga2_f*(8.0*(m1_f+m2_f) + 6.0*m3_f) \
				+ cositot_f*a1_f*CONST_G*m1_f*m2_f*m3_f*pow(8.0*(m1_f+m2_f),-1.0) \
				*(6.0*fme1_f + 18.0*e1_f2*(m1_f*m1_f + m1_f*m2_f + m2_f*m2_f)*cosg1_fd) ) );
	
	Ith(ydot,3) = g1dot_newtonian + f_1PN_1*g1dot_GR_1PN_1 + f_1PN_12*g1dot_GR_1PN_12;


	/* ----------------------------------------------------------------
	* g2dot
	*/ 

	double g2dot_newtonian = 3.0*C2_f*((2.0*cositot_f/Ga1_f)*(2.0 + e1_f2*(3.0 - 5.0*cosg1_fd)) \
			+ (1.0/Ga2_f)*(4.0 + 6.0*e1_f2 + (5.0*cositot_f2 - 3.0)*(2.0 + e1_f2*(3.0 - 5.0*cosg1_fd)))) \
		+ C3_f*e1_f*(sing1_f*sing2_f*(((4.0*e2_f2 + 1.0)/(e2_f*Ga2_f))*10.0*cositot_f*sinitot_f2*(1.0 - e1_f2) \
			- e2_f*((1.0/Ga1_f) + (cositot_f/Ga2_f))*(A_f + 10.0*(3.0*cositot_f2 - 1.0)*(1.0 - e1_f2))) \
			+ cosphi_f*(5.0*B_f*cositot_f*e2_f*((1.0/Ga1_f) + (cositot_f/Ga2_f)) + ((4.0*e2_f2 + 1.0)/(e2_f*Ga2_f))*A_f));
	double g2dot_GR_1PN_2 = 3.0*pow(CONST_C_LIGHT_POW2*a2_f*e2_f2com,-1.0)*pow(CONST_G*m_tot_triple_f/a2_f,3.0/2.0);
	double g2dot_GR_1PN_12 = (CONST_G*pow(CONST_C_LIGHT_POW2,-1.0))*pow( 4.0*pow(a2_f,3.0)*pow(e2_f2com,3.0/2.0)*(m1_f+m2_f),-1.0) \
		*( -2.0*Ga1_f*(8.0*(m1_f+m2_f)+6.0*m3_f)*cositot_f - 3.0*CONST_G*a1_f*m1_f*m2_f*m3_f*pow(8.0*(m1_f+m2_f)*Ga2_f,-1.0) \
				*( fme1_f - 3.0*fme1_f*cositot_f2 + 9.0*e1_f2*(m1_f*m1_f+m1_f*m2_f+m2_f*m2_f)*cosg1_fd*sinitot_f2 ) \
			+ ( (cositot_f/Ga2_f) + (1.0/Ga1_f) )*(Ga1_f*Ga2_f*(8.0*(m1_f+m2_f)+6.0*m3_f) \
				+ CONST_G*a1_f*m1_f*m2_f*m3_f*pow(8.0*(m1_f+m2_f),-1.0)*cositot_f*(6.0*fme1_f \
					+ 18.0*e1_f2*(m1_f*m1_f+m1_f*m2_f+m2_f*m2_f)*cosg1_fd)) );
	Ith(ydot,4) = g2dot_newtonian + f_1PN_2*g2dot_GR_1PN_2 + f_1PN_12*g2dot_GR_1PN_12;		
	

	/* ----------------------------------------------------------------
	* a1dot
	*/ 

	double a1dot_GR_25PN_1 = (-64.0/5.0)*pow(CONST_G,3.0)*pow(CONST_C_LIGHT,-5.0)*m_prod_tot_inner_f*pow(a1_f,-3.0)*pow(e1_f2com,-7.0/2.0)*f_GR_adot1;

	Ith(ydot,5) = f_25PN_1*a1dot_GR_25PN_1;


	/* ----------------------------------------------------------------
	* a2dot
	*/ 

	double a2dot_GR_25PN_2 = (-64.0/5.0)*pow(CONST_G,3.0)*pow(CONST_C_LIGHT,-5.0)*m_tot_inner_f*m3_f*m_tot_triple_f*pow(a2_f,-3.0)*pow(e2_f2com,-7.0/2.0)*f_GR_adot2;

	Ith(ydot,6) = f_25PN_2*a2dot_GR_25PN_2;


	/* ----------------------------------------------------------------
	* cositotdot: due to triple interaction alone!
	*/ 

	double Ga1dot = (-1.0)*Ga1_f*e1_f*(e1dot_newtonian+f_1PN_12*e1dot_GR_1PN_12)/e1_f2com;
	double Ga2dot = (-1.0)*Ga2_f*e2_f*e2dot_newtonian/e2_f2com;
	double cositotdot = (-1.0/(Ga1_f*Ga2_f))*(Ga1dot*(Ga1_f + Ga2_f*cositot_f) + Ga2dot*(Ga2_f + Ga1_f*cositot_f));

	Ith(ydot,7) = cositotdot;


	/* ----------------------------------------------------------------
	* omegaspindot_star_1: not modeled (yet?)
	*/ 

	Ith(ydot,8) = 0.0;


	/* ----------------------------------------------------------------
	* omegaspindot_star_2: not modeled (yet?)
	*/ 	

	Ith(ydot,9) = 0.0;

	return 0;
}

static int froot(realtype t, N_Vector yev, realtype *gout, void *data)
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
