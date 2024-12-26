#ifdef _GSL_
  #include "spmv.h"
#endif
#ifdef _MKL_
  #include "spmv_mkl.h"
#endif
#include <math.h>
#include <stdio.h>

int my_dense(const unsigned int n, const double *restrict mat, double *restrict vec, double *restrict result)
{
  // code your own solver
  unsigned int i, j;
  double tmp = 0.0;
  double *m = (double *)__builtin_assume_aligned(mat, 32);
  double *v = (double *)__builtin_assume_aligned(vec, 32);
  double *r = (double *)__builtin_assume_aligned(result, 32);
#pragma GCC ivdep
  for (i=0; i < n; i++){
  #pragma GCC ivdep
    for (j=0; j < n ; j++){
	tmp += m[i*n+j] * v[j];
    }
	r[i] = tmp;
	tmp=0;
  }

  return 0;
}
