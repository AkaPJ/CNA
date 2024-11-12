#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=5:00
#SBATCH --partition=lpp
OMP_NUM_THREADS=8
mpiexec Ejercicio2hibrido 5 1000000

