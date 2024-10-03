#include "spmv.h"
#include <gsl/gsl_spmatrix.h>
//#include "csr.h"

int my_sparse(const unsigned int n, gsl_spmatrix m, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0, z=0;
  double tmp=0.0;

  for(i=0; i < n; i++){
  	if ((j = m->i[i+1]) != m->i[i]){
		for(j; k < j; k++){
			tmp += m->data[k] * vec[m->p[k]];
		}
		result[z] = tmp;
		tmp = 0;
		++z;
	}
	else {
		result[z] = tmp;
		++z;
	}
  
  }
}
