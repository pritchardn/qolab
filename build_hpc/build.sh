#!/bin/bash

module swap PrgEnv-cray PrgEnv-intel
module unload cray-libsci
module load intel-mkl
module load nlopt
export CRAYPE_LINK_TYPE=dynamic
make
