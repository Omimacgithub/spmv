#include "spmv.h"
#include <gsl/gsl_spmatrix.h>
//#include "csr.h"

int my_sparse(const unsigned int n, /*csr m*/gsl_spmatrix *m, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0, z=0;
  double tmp=0.0;

  for(i=0; i < n; i++){
  	//if ((j = m.row_offsets[i+1]) != m.row_offsets[i]){
  	if ((j = m->p[i+1]) != m->p[i]){
		for(; k < j; k++){
			//tmp += m.values[k] * vec[m.column_indices[k]];
			tmp += m->data[k] * vec[m->i[k]];
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
  return 0;
}
