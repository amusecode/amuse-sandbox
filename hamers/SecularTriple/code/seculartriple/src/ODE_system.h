int fev_triad(realtype t, N_Vector yev, N_Vector ydot, void *data);
int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data);
double f_tides1(double e_p2);
double f_tides2(double e_p2);
double f_tides3(double e_p2);
double f_tides4(double e_p2);
double f_tides5(double e_p2);
double f_25PN_e(double e_p2);
double f_25PN_a(double e_p2);
int froot(realtype t, N_Vector yev, realtype *gout, void *data);

