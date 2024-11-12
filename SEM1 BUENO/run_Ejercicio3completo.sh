#!/bin/bash
#SBATCH--nodes=4
#SBATCH--ntasks=4
#SBATCH--time=5:00
#SBATCH--partition=lpp
mpiexec Ejercicio3completo 8
