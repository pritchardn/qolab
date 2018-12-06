#!/bin/bash

source /opt/intel/parallel_studio_xe_2018/psxevars.sh > /dev/null 2>&1
source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/mklvars.sh intel64
source /opt/intel/compilers_and_libraries_2019.1.144/linux/mpi/intel64/bin/mpivars.sh intel64

#make -C ../../graph_sim_fact
