int fev_triad(realtype t, N_Vector yev, N_Vector ydot, void *data);
int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data);
double f_tides1(double e_p2);
double f_tides2(double e_p2);
double f_tides3(double e_p2);
double f_tides4(double e_p2);
double f_tides5(double e_p2);
double f_25PN_e(double e_p2);
double f_25PN_a(double e_p2);

int froot_delaunay(realtype t, N_Vector yev, realtype *gout, void *data);
double a_out_div_a_in_dynamical_stability(double m1, double m2, double m3, double e_out, double itot);
double a_out_div_a_in_dynamical_stability_mardling_aarseth_01(double m1, double m2, double m3, double e_out, double itot);
double roche_radius_pericenter_eggleton(double rp, double q);
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f);

