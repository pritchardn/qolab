//
// Created by nicholas on 10/12/18.
//
#include <mkl.h>
#include <mathimf.h>
#include "state_evolve.h"


void *initialise_state(MKL_Complex16 *state, machine_spec_t *mach_spec){
    MKL_Complex16 init_value;
    init_value.real = 1/sqrt(mach_spec->space_dimension);
    for(MKL_INT i = 0; i < mach_spec->space_dimension; ++i){
        state[i].real = init_value.real;
    }
}

double evolve(unsigned num_params, const double *x, double grad, qaoa_data_t *meta_spec){
    //Generate new initial state
    MKL_Complex16 *state = mkl_calloc((size_t)meta_spec->machine_spec->space_dimension, sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(state);
    initialise_state(state, meta_spec->machine_spec);
    //for i = 0 to num_params/2
        //diag_expm(uc, v)
        //expm(ub, v)
    //measure (many ways to skin this cat)
    //teardown
    mkl_free(state);
    //Return single value;
    return 0.0;
}
