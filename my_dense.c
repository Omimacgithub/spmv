#ifdef _GSL_
  #include "spmv.h"
#endif
#ifdef _MKL_
  #include "spmv_mkl.h"
#endif
#include <math.h>
#include <stdio.h>

int my_dense(const unsigned int n, const double mat[], double vec[], double result[])
{
  // code your own solver
  int row_offset=n, z=0, i, j, size=n*n;
  double tmp = 0.0;
  
  for (i=0; i < n; i++){
    for (j=0; j < n ; j++){
	tmp += mat[j+i*n] * vec[j];
    }
	result[z] = tmp;
	tmp=0;
	++z;
  }
  //#pragma ivdep
  //#pragma vector always
  /*
  for (int j=0; j < size ; j++){
  //Si haces j%16384 el loop se vectoriza, no se puede, n se conoce en runtime
	//result[j%n] += mat[j] * vec[j];
	
	if((j%n == 0) && (j>0)){
		result[z] = tmp;
		tmp=0;
		++z;
	}
	tmp += mat[j] * vec[j%n];
	
  }
  */
  //result[n-1] = tmp;
  /*
  int j=0;
  while(row_offset <= size){
  	for (; j < row_offset ; j++)
		tmp += mat[j] * vec[j%n];
	result[z] = tmp;
	tmp=0;
	++z;
	row_offset += n;
  }
*/
  
  return 0;
}
