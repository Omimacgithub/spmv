//#include "csr.h"
#include <gsl/gsl_spmatrix.h>

int my_dense(const unsigned int n, const double mat[], double vec[], double result[]);
int my_sparse(const unsigned int n, /*csr m*/gsl_spmatrix *m, double vec[], double result[]);
