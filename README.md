# SpMV: Sparse Matrix-Vector product
## Table of contents:

- [Task](#task)
- [Running the code](#running-the-code)
  - [Makefile](#makefile)
  - [Sbatch job](#sbatch-job)
- [GCC and ICC Benchmark](#gcc-and-icc-benchmark)
  - [Params](#params)
  - [Times](#times)
- [NOTE](#note)

## Task

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

## Running the code

### Makefile
~~~shell
make
~~~

This will create two executables, one named **GSL** (compiled with gcc) and another named **MKL** (with icc).

To only create **GSL**:

~~~shell
make gsl
~~~

To only create **MKL**:

~~~shell
make mkl
~~~

You can use different optimization flags if you want overriding the **OFLAGS** variable:

~~~shell
make mkl OFLAGS="-march=native -fast"
~~~

If you want to remove created object files (.o) **and executable files**:

~~~shell
make cleanall
~~~

### Sbatch job

To perform the benchmark task, I created the file **jub.sbatch** to launch in FT3 the same program several times, each one with a different optimization level.

Just launch the following command:

~~~shell
sbatch jub.sbatch
~~~

## GCC and ICC Benchmark

### Params 

The following experiment was done with the following parameters:
- matrix size: 16384 x 16384
- density (percentage of non-zero elements in the matrix): 10%
- O0 level: -O0
- O2-novec level:
  - GCC: -O2 (no vectorization at this level)
  - ICC: -no-vec -no-simd (-O2 used by default, we need to disable vectorization)
- O3-vec level: -O3 (vectorization enabled) -march=native (let know the compiler what kind of instructions it can emit when generating assembly code)
- Ofast-vec level (gcc): -Ofast -march=native
- fast level (icc): -march=native -fast (if we swap the flags, icc produces the following warning: "overriding '-(null)' with '-march=native'")

### Times

All times are denoted **in miliseconds (ms)**.

| **GSL (GCC)** | O0  | O2-novec | O3-vec | Ofast-vec | Ref             |
| ------------- | --- | -------- | ------ | --------- | --------------- |
| my_dense      | 883 | 587      | 587    | 588       | 148/127/126/130 |
| my_coo        | 141 | 87       | 87     | 87        | 87/87/87/87     |
| my_csr        | 89  | 40       | 40     | 73        | 41/40/40/40     |
| my_csc        | 121 | 35       | 37     | 36        | 37/37/36/36     |

| **MKL (ICC)** | O0  | O2-novec | O3-vec | fast | Ref             |
| ------------- | --- | -------- | ------ | ---- | --------------- |
| my_dense      | 886 | 587      | 588    | 587  | 200/124/128/159 |
| my_coo        | 113 | 87       | 87     | 87   | 117/113/114/113 |
| my_csr        | 89  | 29       | 26     | 25   | 40/38/38/39     |
| my_csc        | 102 | 36       | 34     | 32   | 33/32/33/32     |

## NOTE

In **junkfiles** directory there's stuff no longer needed.
