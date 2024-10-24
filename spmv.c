#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include "timer.h"
#include "spmv.h"
//#include "csr.h"


#define DEFAULT_SIZE 1024
#define DEFAULT_DENSITY 0.25

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
  unsigned int nnz = 0;

  srand(seed);

  for (unsigned int i = 0; i < n * n; i++) {
    if ((rand() % 100) / 100.0 < density) {
      // Get a pseudorandom value between -9.99 e 9.99
      mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
      nnz++;
    } else {
      mat[i] = 0;
    }
  }

  return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
  srand(seed);
  for (unsigned int i = 0; i < size; i++) {
    vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
  }

  return size;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return fabs(x - y) <= epsilon * fabs(x);
  // see Knuth section 4.2.2 pages 217-218
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    if (!is_nearly_equal(ref[i], result[i]))
      return 0;
  }

  return 1;
}

int main(int argc, char *argv[])
{
  int size;        // number of rows and cols (size x size matrix)
  double density;  // aprox. ratio of non-zero values
  //csr m;

  if (argc < 2) {
    size = DEFAULT_SIZE;
    density = DEFAULT_DENSITY;
  } else if (argc < 3) {
    size = atoi(argv[1]);
    density = DEFAULT_DENSITY;
  } else {
    size = atoi(argv[1]);
    density = (double) atoi(argv[2]) / 100.0;
  }

  gsl_spmatrix *m = gsl_spmatrix_compress(gsl_spmatrix_alloc(size, size), GSL_SPMATRIX_CSR); //gsl_spmatrix in CSR format
  gsl_spmatrix *src = gsl_spmatrix_alloc(size, size);	   // gsl_spmatrix 
  double *mat, *vec, *refsol, *mysol;

  mat = (double *) malloc(size * size * sizeof(double));
  vec = (double *) malloc(size * sizeof(double));
  refsol = (double *) malloc(size * sizeof(double));
  mysol = (double *) malloc(size * sizeof(double));

  unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
  populate_vector(vec, size, 2);

  printf("Matriz size: %d x %d (%d elements)\n", size, size, size*size);
  printf("%d non-zero elements (%.2lf%%)\n\n", nnz, (double) nnz / (size*size) * 100.0);

  //
  // Dense computation using CBLAS (eg. GSL's CBLAS implementation)
  //
  printf("Dense computation\n----------------\n");

  timeinfo start, now;
  timestamp(&start);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

  timestamp(&now);

  printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));

  //
  // Using your own dense implementation
  //
  timestamp(&start);

  my_dense(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  //
  // Let's try now SpMV: Sparse Matrix - Dense Vector computation
  //

  // Convert mat to a sparse format: CSR
  // Use the gsl_spmatrix struct as datatype
  //m = intoCSR(size, mat);
  double value;
  int i, j;
  for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            value = mat[i * size + j];
            if (value != 0.0) {
                gsl_spmatrix_set(src, i, j, value);
            }
        }
    }

  //Store CSR matrix in m
  gsl_spmatrix_csr(m, src);
  //
  // csr computation using GSL's sparse algebra functions
  //
  
  gsl_vector *x = gsl_vector_alloc(size);

  for (i = 0; i < size; i++) {
  	gsl_vector_set(x, i, vec[i]);
  }

  gsl_vector *y = gsl_vector_alloc(size);
  gsl_vector_set_zero(y);

  //Matrix-vector operation in the form: y = alpha*m*x + beta*y
  //alpha = 1 and beta = 0

  timestamp(&start);

  //Result stored in y gsl_vector
  gsl_spblas_dgemv(CblasNoTrans, 1.0, m, x, 0.0, y);

  timestamp(&now);
  printf("Time taken by GSL (CSR) sparse matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  for(i=0; i < size; i++){
	mysol[i] = gsl_vector_get(y, i);
  }

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  //
  // Your own csr implementation
  //

  // Compare times (and computation correctness!)
  timestamp(&start);

  my_csr(size, m, vec, mysol);

  timestamp(&now);
  printf("Time taken by my csr matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

//Free and reinitialize stuff

  gsl_spmatrix_free(m);

  m = gsl_spmatrix_compress(gsl_spmatrix_alloc(size, size), GSL_SPMATRIX_CSC); //gsl_spmatrix in CSC format

  //Store CSC matrix in m
  gsl_spmatrix_csc(m, src);
  //
  // csc computation using GSL's sparse algebra functions
  //
  
  //Matrix-vector operation in the form: y = alpha*m*x + beta*y
  //alpha = 1 and beta = 0

  timestamp(&start);

  //Result stored in y gsl_vector
  gsl_spblas_dgemv(CblasNoTrans, 1.0, m, x, 0.0, y);

  timestamp(&now);
  printf("Time taken by GSL (CSC) sparse matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  for(i=0; i < size; i++){
	mysol[i] = gsl_vector_get(y, i);
  }

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  //
  // Your own csc implementation
  //

  // Compare times (and computation correctness!)
  timestamp(&start);

  my_csc(size, m, vec, mysol);

  timestamp(&now);
  printf("Time taken by my csc matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  //
  // COO computation using GSL's sparse algebra functions
  //
  
  //Matrix-vector operation in the form: y = alpha*m*x + beta*y
  //alpha = 1 and beta = 0

  timestamp(&start);

  //Result stored in y gsl_vector
  gsl_spblas_dgemv(CblasNoTrans, 1.0, src, x, 0.0, y);

  timestamp(&now);
  printf("Time taken by GSL (COO) sparse matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  for(i=0; i < size; i++){
	mysol[i] = gsl_vector_get(y, i);
  }

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  //
  // Your own coo implementation
  //

  // Compare times (and computation correctness!)
  timestamp(&start);

  my_coo(size, src, vec, mysol);

  timestamp(&now);
  printf("Time taken by my coo matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  // Free resources
  free(mat);
  free(vec);
  free(refsol);
  free(mysol);
  gsl_spmatrix_free(m);
  gsl_spmatrix_free(src);
  gsl_vector_free(x);
  gsl_vector_free(y);
  return 0;
}
