#include "spmv.h"
#include <gsl/gsl_spmatrix.h>

int my_csc(const unsigned int n, gsl_spmatrix *m, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0;
  for (i=0; i < n; i++)
  	result[i]=0;
  // p = column_offsets
  // i = row_indices
  // m->data = values
  // Docs: https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix_csr
  // n siendo el nº de columnas de la matriz
  for(i=0; i < n; i++){
  	if ((j = m->p[i+1]) != m->p[i]){
		for(; k < j; k++){
			//tmp += m.values[k] * vec[m.column_indices[k]];
			//Multiplica los valores no nulos de la matriz (en las posiciones de column_indices) con el vector
			result[m->i[k]] += m->data[k] * vec[i];
		}
	}
  }
  return 0;
}
