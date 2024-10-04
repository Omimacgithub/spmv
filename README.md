# SpMV: Sparse Matrix-Vector product

Use this code skeleton for the associated tasks in HPCTools.

This code is based on the use of GSL (GNU Scientific Library) for the
implementation of the baseline operations used for comparison:
- dense matrix-vector product: `cblas_dgemv()`, you need to link against *libgslcblas*
- sparse matrix-vector product: `gsl_spblas_dgemv()`, you need to link against *libgsl*

The dense product, cblas_dgemv(), can be found in other CBLAS
implementation. You just need to change the library to be linked,
eg. `-lopenblas` instead of `-lgslcblas`

The basetype in GSL for working with sparse matrices is `gsl_spmatrix`.
GSL also provides functions to help convert you dense matrices into a sparse format.

This code is intended to use in the course HPCTools

## Runing the code

~~~shell
make
~~~

This will create two executables, one named **release** (with -O2 CFLAG) and another named **debug** (with -O0 and -g3 CFLAGs).

To only create **release**:

~~~shell
make release
~~~

To only create **debug**:

~~~shell
make debug
~~~

If you want to remove created object files (.o):

~~~shell
make clean
~~~

If you want to remove created object files (.o) **and executable files**:

~~~shell
make cleanall
~~~

In **junkfiles** directory there's stuff no longer needed.

