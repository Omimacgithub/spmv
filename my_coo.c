#ifdef _GSL_
#include "spmv.h"
#include <stddef.h>
#include <gsl/gsl_spmatrix.h>
#endif
#ifdef _MKL_
#include <mkl.h>
#include "spmv_mkl.h"
#endif

#ifdef _GSL_
int my_coo(const unsigned int n, gsl_spmatrix *m, double vec[], double result[])

{
  // code your own solver
  int i, size = m->nz;
  //m->data = values
  //m->p = column_indices
  //m->i = row_indices
  //Docs: https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix_csr
  for (i=0; i < n; i++)
  	result[i]=0;
  for(i=0; i < size; i++)
	result[m->i[i]]	+= m->data[i] * vec[m->p[i]];
  
  return 0;
}
#endif

#ifdef _MKL_
int my_coo(const unsigned int n, const unsigned int nnz, MKL_INT *rows_indx, MKL_INT *cols_indx, const double *values, double vec[], double result[])
{
  // code your own solver
  //sparse_index_base_t indexing;
  //MKL_INT nrows, ncols, nnz;
  //const MKL_INT *rows_indx, *cols_indx;
  //const double *values;

  
  //You CANNOT access to the elements of sparse_matrix_t by hand (are opaque)
  //mkl_sparse_d_export_coo(m, &indexing, &nrows, &ncols, &nnz, &rows_indx, &cols_indx, &values);

  int i;
  for(i=0; i < n; i++)
  	result[i]=0;

  for(i=0; i < nnz; i++)
	result[rows_indx[i]] += values[i] * vec[cols_indx[i]];
  
  return 0;
}
#endif
