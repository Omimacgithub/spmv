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
#endif
#ifdef _MKL_
//int my_csc(const unsigned int n, sparse_matrix_t *m, double vec[], double result[])
int my_csc(const unsigned int n, MKL_INT *cols_start, MKL_INT *rows_indx, const double *values, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0;
  for (i=0; i < n; i++)
  	result[i]=0;
  /*
   * sparse_index_base_t indexing;
  MKL_INT nrows, ncols, nnz; 
  MKL_INT *cols_start, *cols_end, *rows_indx;
  double *values;

  mkl_sparse_d_export_csc(*m, &indexing, &nrows, &ncols, &cols_start, &cols_end, &rows_indx, &values);
  */
  for(i=0; i < n; i++){
  	if ((j = cols_start[i+1]) != cols_start[i]){
		for(; k < j; k++){
			result[rows_indx[k]] += values[k] * vec[i];
		}
	}
  }
  return 0;
}
#endif
