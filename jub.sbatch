#!/bin/bash -l
#Docs: https://cesga-docs.gitlab.io/ft3-user-guide/batch_examples.html
#SBATCH -n 1                     #Nº of MPI/process tasks (1)
#SBATCH --time=0-00:10:00
#SBATCH --mem=8G
#SBATCH --exclusive              #Allocate an exclusive node (for execution time fidelity)

SIZE=16384
DENSITY=10

#Load icc and mkl
module load cesga/2022 iimkl

make cleanall
make OFLAGS="-O0"

echo "-------------------------------"
echo "--------Results with -O0-------"
echo "-------------------------------"

echo "-------------GSL---------------"
./gsl $SIZE $DENSITY

echo "-------------MKL---------------"
./mkl $SIZE $DENSITY

make cleanall
make gsl OFLAGS="-O2"
make mkl OFLAGS="-no-vec -no-simd"
echo "-------------------------------"
echo "-----Results with -O2-novec----"
echo "-------------------------------"

echo "-------------GSL---------------"
./gsl $SIZE $DENSITY

echo "-------------MKL---------------"
./mkl $SIZE $DENSITY

make cleanall
make OFLAGS="-O3 -march=native"
echo "-------------------------------"
echo "--------Results with -O3-------"
echo "-------------------------------"

echo "-------------GSL---------------"
./gsl $SIZE $DENSITY

echo "-------------MKL---------------"
./mkl $SIZE $DENSITY

make cleanall
make gsl OFLAGS="-Ofast -march=native"
make mkl OFLAGS="-march=native -fast"

echo "-------------------------------"
echo "------Results with -Ofast------"
echo "-------------------------------"

echo "-------------GSL---------------"
./gsl $SIZE $DENSITY

echo "-------------MKL---------------"
./mkl $SIZE $DENSITY
