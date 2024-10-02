#include "spmv.h"
#include "csr.h"

int my_sparse(const unsigned int n, csr m, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0, z=0;
  double tmp=0.0;

  for(i=0; i < n; i++){
  	if ((j = m.row_offsets[i+1]) != m.row_offsets[i]){
		for(j; k < j; k++){
			tmp += m.values[k] * vec[m.column_indices[k]];
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
