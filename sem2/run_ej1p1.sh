#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=5:00
#SBATCH --partition=lpp
OMP_NUM_THREADS=4 ./ej1programa1 500
