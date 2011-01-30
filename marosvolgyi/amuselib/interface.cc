#include "worker_code.h"
#ifdef __cplusplus
extern "C" {
#endif
int asl_add(double x, double y, double *z);
#ifdef __cplusplus
}
#endif

int add(double x, double y, double *z)
{
  asl_add(x, y, z);
  return 0;
}
