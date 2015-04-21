#include "main_code.h"


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
