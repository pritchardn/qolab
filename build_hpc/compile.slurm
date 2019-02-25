#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:10:00
#SBATCH --account=pawsey0309
#SBATCH --export=NONE
 
module swap PrgEnv-cray PrgEnv-intel
module load intel-mkl
module load nlopt
 
srun --export=all cc -O2 -xCORE-AVX2 ../src/main.c ../src/qaoa.c ../src/ub.c ../src/globals.c ../src/uc.c ../src/problem_code.c ../src/state_evolve.c ../src/reporting.c ../src/matrix_expm.c ../src/graph_utils.c ../src/qaoa.h ../src/ub.h ../src/globals.h ../src/uc.h ../src/problem_code.h ../src/state_evolve.h ../src/reporting.h ../src/matrix_expm.h ../src/graph_utils.h
