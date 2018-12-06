#include "qaoa.h"
#include "ub.h"
#include <stdio.h>


void qaoa(machine_spec_t *mach_spec, optimisation_spec_t *opt_spec, run_spec_t *run_spec){
    qaoa_data_t meta_spec;
    meta_spec.machine_spec = mach_spec;
    meta_spec.run_spec = run_spec;

    printf("%d\n", meta_spec.machine_spec->num_qubits);

    sparse_matrix_t ub;
    generate_ub(mach_spec->num_qubits, &ub);
    mkl_sparse_print(&ub, stdout);
    mkl_sparse_destroy(ub);

    //Initialise UC
    //Initialise UB
    //Trial feval
    //Teardown
}
