#include <stddef.h>
#ifdef _GSL_
#include <gsl/gsl_spmatrix.h>
#include "spmv.h"
#endif
#ifdef _MKL_
#include <mkl.h>
#include "spmv_mkl.h"
#endif
//#include "csr.h"

#ifdef _GSL_
int my_csr(const unsigned int n, /*csr m*/gsl_spmatrix *m, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0, z=0;
  double tmp=0.0;
  //m->data = values
  //m->p = row_offsets
  //m->i = column_indices
  //Docs: https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix_csr

  for(i=0; i < n; i++){
  	//if ((j = m.row_offsets[i+1]) != m.row_offsets[i]){
  	if ((j = m->p[i+1]) != m->p[i]){
		for(; k < j; k++){
			//tmp += m.values[k] * vec[m.column_indices[k]];
			//Multiplica los valores no nulos de la matriz (en las posiciones de column_indices) con el vector
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
#endif
#ifdef _MKL_
int my_csr(const unsigned int n, sparse_matrix_t *m, double vec[], double result[])
{
  // code your own solver
  int i, j=0, k=0, z=0;
  double tmp=0.0;
  sparse_index_base_t indexing;
  MKL_INT nrows, ncols, nnz; 
  MKL_INT *rows_start, *rows_end, *cols_indx;
  double *values;

  //You CANNOT access to the elements of sparse_matrix_t by hand (are opaque)
  mkl_sparse_d_export_csr(*m, &indexing, &nrows, &ncols, &rows_start, &rows_end, &cols_indx, &values);


  for(i=0; i < n; i++){
  	if ((j = rows_start[i+1]) != rows_start[i]){
		for(; k < j; k++){
			tmp += values[k] * vec[cols_indx[k]];
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
#endif
