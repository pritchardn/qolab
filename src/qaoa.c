#include "qaoa.h"
#include "ub.h"
#include <stdio.h>


void qaoa_teardown(qaoa_data_t *meta_spec){
    mkl_sparse_destroy(meta_spec->ub);
}

void qaoa(machine_spec_t *mach_spec, optimisation_spec_t *opt_spec, run_spec_t *run_spec){
    qaoa_data_t meta_spec;
    meta_spec.machine_spec = mach_spec;
    meta_spec.run_spec = run_spec;

    printf("%d\n", meta_spec.machine_spec->num_qubits);

    //Initialise UC

    //Initialise UB
    generate_ub(&meta_spec);
    //mkl_sparse_print(&meta_spec.ub, stdout);
    //Trial feval
    //Teardown
    qaoa_teardown(&meta_spec);
}
