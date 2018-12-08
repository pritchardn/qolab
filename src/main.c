//
// Created by nicholas on 4/12/18.
//

#include "main.h"
#include "qaoa.h"
#include <math.h>

int main(int argc, char *argv[]){

    run_spec_t run_spec;
    run_spec.correct = true;
    run_spec.report = false;
    run_spec.timing = true;

    optimisation_spec_t opt_spec;
    opt_spec.ftol = 1e-8;
    opt_spec.xtol = 1e-8;
    opt_spec.nlopt_method = 28;

    machine_spec_t mach_spec;
    mach_spec.num_qubits = 20;
    mach_spec.P = 1;
    mach_spec.space_dimension = (MKL_INT)pow(2, mach_spec.num_qubits);

    qaoa(&mach_spec, &opt_spec, &run_spec);

    return 0;
}
