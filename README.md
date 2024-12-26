# SpMV: Sparse Matrix-Vector product
## Table of contents:

- [Authors](#authors)
- [Task](#task)
- [Running the code](#running-the-code)
  - [Makefile](#makefile)
  - [Sbatch job](#sbatch-job)
- [GCC and ICC Benchmark](#gcc-and-icc-benchmark)
  - [Params](#params)
  - [Times](#times)
- [NOTE](#note)

## Authors

- Omar Montenegro Macía
- Emilio J. Padrón González (provider of the [code skeleton](https://gitlab.citic.udc.es/emilio.padron/spmv))

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

You can use different optimization flags if you want by overriding the **OFLAGS** variable:

~~~shell
make mkl OFLAGS="-march=native -fast"
~~~

If you want to remove created object files (.o) **and executable files**:

~~~shell
make cleanall
~~~

### Sbatch job

To perform the benchmark task, I created the file **jub.sbatch** to launch in FT3 the same program several times, each one with a different optimization level. The script launch the experiment **3 times** (the times the experiment is executed is controlled by the variable **num_executions**) and calculate the mean of all execution times.

Just launch the following command:

~~~shell
sbatch jub.sbatch
~~~

## GCC and ICC Benchmark

### Params 

The experiment was done with the following parameters:
- matrix size: 16384 x 16384
- density (percentage of non-zero elements in the matrix): 10%
- O0 level: -O0
- O2-novec level:
  - GCC: -O2 (no vectorization at this level)
  - ICC: -no-vec -no-simd (-O2 used by default, we need to disable vectorization)
- O3-vec level: -O3 (vectorization enabled) -march=native (let know the compiler what kind of instructions it can emit when generating assembly code)
- Ofast-vec level (gcc): -Ofast -march=native
- fast level (icc): -march=native -fast (if we swap the flags order, icc produces the following warning: "overriding '-(null)' with '-march=native'")
- Libraries: -lm
  - GSL: -lgslcblas (for dense matrix vector multiplication) -lgsl (for sparse matrix vector multiplication)
  - MKL: -lmkl
- gcc version: 11.2.1 20220115 (Gentoo 11.2.1_p20220115 p4)
- icc version: 2021.10.0 Build 20230609_000000

### Times

The times reflected in the 2 following tables are the **mean of executing the same program 3 times**. All times are denoted **in miliseconds (ms)**.

The **Ref** column shows the execution time of GSL or MKL library function (already compiled and optimized).

| **GSL (GCC)** | O0  | O2-novec | O3-vec | Ofast-vec | Ref             |
| ------------- | --- | -------- | ------ | --------- | --------------- |
| my_dense      | 879.33 | 587      | 587    | 587       | 398.66 |
| my_coo        | 140.33 | 86.33       | 87     | 87        | 87     |
| my_csr        | 88.66  | 40       | 40     | 73.66        | 40.33     |
| my_csc        | 119.66 | 35.66       | 36     | 35.66        | 36     |

| **MKL (ICC)** | O0  | O2-novec | O3-vec | fast | Ref             |
| ------------- | --- | -------- | ------ | ---- | --------------- |
| my_dense      | 888 | 588      | 587.33    | 590.33  | 189.66 |
| my_coo        | 112.33 | 87       | 86.66     | 86.66   | 114.33 |
| my_csr        | 88.33  | 29       | 25.66     | 25   | 39     |
| my_csc        | 102.66 | 35.66       | 34.66     | 32   | 32  |

## NOTE

In **junkfiles** directory there's stuff no longer needed.
