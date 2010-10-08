#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>



int main(int argc, char *argv[]) {
  mpfr_t x, y;

  mpfr_inits2(100, x, y, (mpfr_ptr)0);
  mpfr_set_d(x, 10.0, GMP_RNDN);
  mpfr_printf ("variable x with %.50Rf \n", x);
  return 0;
}
