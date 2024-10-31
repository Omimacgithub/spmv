#!/bin/bash -l
#Docs: https://cesga-docs.gitlab.io/ft3-user-guide/batch_examples.html
#SBATCH -N 1                      #Nº de nodos (1)
#SBATCH -n 1                     #Nº de tareas MPI/procesos (1)
#SBATCH -c 1                      #Nº de cores por proceso/tarea (1)
#SBATCH --ntasks-per-node=1	  #Nº de tareas por nodo (1)
#SBATCH --time=0-00:10:00
#SBATCH --mem=200G

module load cesga/2022 iimkl

./mkl 16384
