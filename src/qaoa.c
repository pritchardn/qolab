#include "qaoa.h"
#include "ub.h"
#include "uc.h"
#include "problem_code.h"
#include <stdio.h>


void qaoa_teardown(qaoa_data_t *meta_spec){
    mkl_sparse_destroy(meta_spec->ub);
    mkl_free(meta_spec->uc);
}

void qaoa(machine_spec_t *mach_spec, cost_data_t *cost_data, optimisation_spec_t *opt_spec, run_spec_t *run_spec){
    qaoa_data_t meta_spec;
    qaoa_statistics_t statistics;
    meta_spec.qaoa_statistics = &statistics;
    meta_spec.machine_spec = mach_spec;
    meta_spec.run_spec = run_spec;

    printf("%d\n", meta_spec.machine_spec->num_qubits);

    //Initialise UC
    meta_spec.uc = mkl_malloc((size_t)pow(2, meta_spec.machine_spec->num_qubits) * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    generate_uc(&meta_spec, cost_data, Cx);
    //Initialise UB
    generate_ub(&meta_spec);
    //mkl_sparse_print(&meta_spec.ub, stdout);
    //Trial feval
    //Teardown
    qaoa_teardown(&meta_spec);
}
