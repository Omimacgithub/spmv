#ifdef _GSL_
#include "spmv.h"
#include <gsl/gsl_spmatrix.h>
#endif
#ifdef _MKL_
#include <mkl.h>
#include "spmv_mkl.h"
#endif

#ifdef _GSL_
int my_csc(const unsigned int n, gsl_spmatrix *m, double vec[], double result[])
{
  // code your own solver
  unsigned int i, j=0, k=0;
  double *Md = m->data;
  int *Mp = m->p, *Mi = m->i;

  for (i=0; i < n; i++)
  	result[i]=0;
  // p = column_offsets
  // i = row_indices
  // m->data = values
  // Docs: https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix_csr
  // n siendo el nº de columnas de la matriz
  for(i=0; i < n; i++){
  	j = Mp[i+1];
	for(k = Mp[i]; k < j; k++){
		//Multiplica los valores no nulos de la matriz (en las posiciones de column_indices) con el vector
		result[Mi[k]] += Md[k] * vec[i];
		}
	}

  return 0;
}
#endif
#ifdef _MKL_
int my_csc(const unsigned int n, MKL_INT *cols_start, MKL_INT *rows_indx, const double *values, double vec[], double result[])
{
  // code your own solver
  unsigned int i, j=0, k=0;
  for (i=0; i < n; i++)
  	result[i]=0;
  for(i=0; i < n; i++){
  	j = cols_start[i+1];
		for(k=cols_start[i] ; k < j; k++){
			result[rows_indx[k]] += values[k] * vec[i];
		}
	}
  
  return 0;
}
#endif
