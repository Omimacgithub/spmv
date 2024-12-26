#!/bin/bash -l
#Docs: https://cesga-docs.gitlab.io/ft3-user-guide/batch_examples.html
#SBATCH -n 1                     #Nº of MPI/process tasks (1)
#SBATCH -o mkl_%j.o
#SBATCH --time=0-00:20:00
#SBATCH --mem=8G
#SBATCH --exclusive              #Allocate an exclusive node (for execution time fidelity)

SIZE=16384
DENSITY=10

#Load icx and mkl
module load cesga/2020 intel imkl
icx --version

# Number of executions
num_executions=3

# Declare an associative array to store total times for each metric
declare -A execution_times_O0
declare -A execution_times_O2
declare -A execution_times_O3
declare -A execution_times_fast

# Initialize the array for each metric
metrics=("CBLAS (MKL) dense computation" "my dense matrix-vector product (MKL)" \
         "MKL (CSR) sparse matrix - vector product" "my csr matrix (MKL) - vector product" \
         "MKL (CSC) sparse matrix - vector product" "my csc matrix (MKL) - vector product" \
         "MKL (COO) sparse matrix - vector product" "my coo matrix (MKL) - vector product")
# Initialize counts to zero
for metric in "${metrics[@]}"; do
    execution_times_O0["$metric"]=0
    execution_times_O2["$metric"]=0
    execution_times_O3["$metric"]=0
    execution_times_fast["$metric"]=0
done


# Execute the C program multiple times
for (( i=1; i<=num_executions; i++ )); do
    make cleanall
    make mkl OFLAGS="-O0" IVEC=""

    output=$(./mkl $SIZE $DENSITY)

   echo "-------------------------------"
   echo "--------Results with -O0-------"
   echo "-------------------------------"
   echo $output

    # Read output line by line
    while IFS= read -r line; do
        # Use regex to capture the metric and time
        if [[ $line =~ Time\ taken\ by\ (.*?):\ ([0-9]+)\ ms ]]; then
            metric="${BASH_REMATCH[1]}"
            time="${BASH_REMATCH[2]}"
            execution_times_O0["$metric"]=$((execution_times_O0["$metric"] + time))
        fi
    done <<< "$output"

    make cleanall
    make mkl OFLAGS="-no-vec -no-simd" IVEC=""

    output=$(./mkl $SIZE $DENSITY)

   echo "-------------------------------"
   echo "-----Results with -O2-novec----"
   echo "-------------------------------"
   echo $output

    while IFS= read -r line; do
        # Use regex to capture the metric and time
        if [[ $line =~ Time\ taken\ by\ (.*?):\ ([0-9]+)\ ms ]]; then
            metric="${BASH_REMATCH[1]}"
            time="${BASH_REMATCH[2]}"
            execution_times_O2["$metric"]=$((execution_times_O2["$metric"] + time))
        fi
    done <<< "$output"

    make cleanall
    make mkl OFLAGS="-O3 -march=native" 

    output=$(./mkl $SIZE $DENSITY)

   echo "-------------------------------"
   echo "--------Results with -O3-------"
   echo "-------------------------------"
   echo $output

    while IFS= read -r line; do
        # Use regex to capture the metric and time
        if [[ $line =~ Time\ taken\ by\ (.*?):\ ([0-9]+)\ ms ]]; then
            metric="${BASH_REMATCH[1]}"
            time="${BASH_REMATCH[2]}"
            execution_times_O3["$metric"]=$((execution_times_O3["$metric"] + time))
        fi
    done <<< "$output"

    make cleanall
    make mkl OFLAGS="-march=native -fast"

    output=$(./mkl $SIZE $DENSITY)

   echo "-------------------------------"
   echo "-------Results with -Ofast-----"
   echo "-------------------------------"
   echo $output
   
    while IFS= read -r line; do
        # Use regex to capture the metric and time
        if [[ $line =~ Time\ taken\ by\ (.*?):\ ([0-9]+)\ ms ]]; then
            metric="${BASH_REMATCH[1]}"
            time="${BASH_REMATCH[2]}"
            execution_times_fast["$metric"]=$((execution_times_fast["$metric"] + time))
        fi
    done <<< "$output"
done

# Calculate and print mean times
for metric in "${metrics[@]}"; do
	echo "------------------------"
	echo "-----------O0-----------"
	echo "------------------------"
        mean_time=$(echo "scale=2; ${execution_times_O0["$metric"]} / $num_executions" | bc)
        echo "Mean time for $metric: $mean_time ms"

	echo "------------------------"
	echo "-----------O2-----------"
	echo "------------------------"
        mean_time=$(echo "scale=2; ${execution_times_O2["$metric"]} / $num_executions " | bc)
        echo "Mean time for $metric: $mean_time ms"

	echo "------------------------"
	echo "-----------O3-----------"
	echo "------------------------"
        mean_time=$(echo "scale=2; ${execution_times_O3["$metric"]} / $num_executions" | bc)
        echo "Mean time for $metric: $mean_time ms"

	echo "------------------------"
	echo "---------fast-----------"
	echo "------------------------"

        mean_time=$(echo "scale=2; ${execution_times_fast["$metric"]} / $num_executions " | bc)
        echo "Mean time for $metric: $mean_time ms"
done

