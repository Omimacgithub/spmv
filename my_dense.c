#ifdef _GSL_
  #include "spmv.h"
#endif
#ifdef _MKL_
  #include "spmv_mkl.h"
#endif
#include <math.h>
#include <stdio.h>

int my_dense(const unsigned int n, const double  mat[], double  vec[], double  result[])
{
  // code your own solver
  int row_offset=n, z=0, i, j, size=n*n;
  double tmp = 0.0;
  
  #pragma GCC ivdep
 // #pragma omp simd reduction(+:tmp)
  for (i=0; i < n; i++){
  #pragma GCC ivdep
    for (j=0; j < n ; j++){
	tmp += mat[j+i*n] * vec[j];
    }
	result[z] = tmp;
	tmp=0;
	++z;
  }
  
  return 0;
}
