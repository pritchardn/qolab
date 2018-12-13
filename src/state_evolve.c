//
// Created by nicholas on 10/12/18.
//
#include <mkl.h>
#include <mathimf.h>
#include "state_evolve.h"
#include "matrix_expm.h"


void extract_hamiltonian_double(MKL_Complex16 *uc, double *hamiltonian, MKL_INT space_dimension){
    for (int i = 0; i < space_dimension; ++i) {
        hamiltonian[i] = -uc[i].imag;
    }
}

void *initialise_state(MKL_Complex16 *state, machine_spec_t *mach_spec){
    MKL_Complex16 init_value;
    init_value.real = 1/sqrt(mach_spec->space_dimension);
    for(MKL_INT i = 0; i < mach_spec->space_dimension; ++i){
        state[i].real = init_value.real;
    }
}

void compute_probabilities(MKL_Complex16 * state, qaoa_data_t *meta_spec, double *output){
    MKL_Complex16 *z_probabilities = mkl_malloc(meta_spec->machine_spec->space_dimension * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(z_probabilities);
    //double result = 0.0;
    //cblas_zdotc_sub(space_dimension, state, 1, state, 1, &result);
    //printf("%f\n", result);
    vzMulByConj(meta_spec->machine_spec->space_dimension, state, state, z_probabilities);
    vzAbs(meta_spec->machine_spec->space_dimension, z_probabilities, output);
    mkl_free(z_probabilities);
}

double expectation_value(MKL_Complex16 *state, qaoa_data_t *meta_spec){
    double *hamiltonian = NULL;
    double expectation;
    MKL_INT space_dimension = meta_spec->machine_spec->space_dimension;

    double *probabilities = mkl_malloc(space_dimension * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    compute_probabilities(state, meta_spec, probabilities);

    hamiltonian = mkl_malloc(space_dimension * sizeof(double), DEF_ALIGNMENT);
    extract_hamiltonian_double(meta_spec->uc, hamiltonian, space_dimension);
    expectation = cblas_ddot(space_dimension, probabilities, 1, hamiltonian, 1);
    /*int ind_max = 0;
    double max_prob = (double) -INFINITY;
    for (int i = 0; i < space_dimension; ++i) {
        if (probabilities[i] > max_prob) {
            max_prob = probabilities[i];
            //ind_max = i;
        }
    }*/
    mkl_free(hamiltonian);
    return expectation;
}

double evolve(unsigned num_params, const double *x, double *grad, qaoa_data_t *meta_spec){
    double result = 0.0;
    int P = meta_spec->machine_spec->P;
    //Generate new initial state
    MKL_Complex16 *state = mkl_calloc((size_t)meta_spec->machine_spec->space_dimension, sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(state);
    initialise_state(state, meta_spec->machine_spec);
    for(int i = 0; i < num_params/2; ++i){
        spmatrix_expm_z_diag(meta_spec->uc, x[i], meta_spec->machine_spec->space_dimension, state);
        spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16){x[i+P], 0.0}, (MKL_Complex16){0.0, -meta_spec->machine_spec->num_qubits},(MKL_Complex16){0.0, meta_spec->machine_spec->num_qubits}, meta_spec->machine_spec->space_dimension);
    }
    //measure (many ways to skin this cat)
    if(meta_spec->run_spec->sampling){
        //TODO: Support compressed sampling
    } else {
        result = expectation_value(state, meta_spec);
    }
    //teardown
    mkl_free(state);
    //Return single value;
    return result;
}
