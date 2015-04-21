int compute_triad_vectors_from_orbital_elements(double m1, double m2, double m3, double a1, double a2, double e1, double e2, double INCL1, double INCL2, double AP1, double AP2, double LAN1, double LAN2,
    double e1_vec[3], double e2_vec[3], double h1_vec[3], double h2_vec[3], double q1_vec[3], double q2_vec[3]);
int compute_orbital_elements_from_triad_vectors(double m1, double m2, double m3,
    double e1_vec[3],double e2_vec[3],double h1_vec[3],double h2_vec[3],double q1_vec[3],double q2_vec[3],
    double *a1_out, double *a2_out, double *e1_out, double *e2_out, double *INCL1_out, double *INCL2_out, double *AP1_out, double *AP2_out, double *LAN1_out, double *LAN2_out);
    
void cross3(double a[3], double b[3], double result[3]);
double norm3(double v[3]);
double norm3_squared(double v[3]);
double dot3(double a[3], double b[3]);
