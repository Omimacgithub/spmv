#include "spmv.h"
#include <math.h>

int my_dense(const unsigned int n, const double mat[], double vec[], double result[])
{
  // code your own solver
  int rows=0, j=0, z=0;
  double tmp = 0.0;
  //Assuming an NxN matrix 
  rows = (int) sqrt(n);
  while(rows <= n ){
  	for (j=j; j < rows ; ++j){
		tmp = tmp + mat[j] * vec[j];

	}
	result[z] = tmp;
	tmp=0;
	++z;
	rows += j;
  }
}
