#include "spmv.h"
#include <gsl/gsl_spmatrix.h>
//#include "csr.h"

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
