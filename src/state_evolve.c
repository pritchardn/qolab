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
    //TODO: Unit test for normalisation
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
    //A potential measurement method, finds the most likely outcome.
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

MKL_INT binary_search(const double *array, const double target, MKL_INT size) {
    MKL_INT first = 0;
    MKL_INT last = size - 1;
    MKL_INT middle;
    while (first <= last) {
        middle = ((first + last) / 2);
        if (array[middle] < target) {
            first = middle + 1;
        } else if (array[middle] > target) {
            last = middle - 1;
        } else {
            return middle;
        }
    }
    if (last < 0) {
        return 0;
    } else if (first > size - 1) {
        return size - 1;
    } else {
        return (first < last ? first : last + 1);
    }
}

void weighted_choice(const MKL_INT *vals, const double *weights, MKL_INT size, run_spec_t *run_spec, double *result) {
    MKL_INT temp;
    double target;
    for (int i = 0; i < run_spec->num_samples; ++i) {
        target = rand() / (double) RAND_MAX;
        temp = binary_search(weights, target, size);
        result[i] = (vals[temp]);
    }
}

MKL_INT build_interval_array(qaoa_data_t *meta_spec, const double *probabilities, MKL_INT *result_vals,
                             double *result_weights) {
    double *hamiltonian;
    double *sum_vals;
    bool *set_flag;
    MKL_INT space_dimension = meta_spec->machine_spec->space_dimension;
    MKL_INT nnz = 0;
    MKL_INT num_vals = meta_spec->cost_data->cx_range;
    hamiltonian = mkl_malloc(space_dimension * sizeof(double), DEF_ALIGNMENT);
    check_alloc(hamiltonian);
    sum_vals = mkl_calloc((size_t) num_vals, sizeof(double), DEF_ALIGNMENT);
    check_alloc(sum_vals);
    set_flag = mkl_calloc((size_t) num_vals, sizeof(bool), DEF_ALIGNMENT);
    check_alloc(set_flag);

    extract_hamiltonian_double(meta_spec->uc, hamiltonian, space_dimension);
    for (MKL_INT i = 0; i < space_dimension; ++i) {
        if (!set_flag[((int) hamiltonian[i]) - 1]) {
            result_vals[nnz] = ((int) hamiltonian[i]);
            nnz++;
            set_flag[(int) hamiltonian[i] - 1] = true;
        }
        sum_vals[(int) hamiltonian[i] - 1] += probabilities[i];
    }
    vdPackV(nnz, sum_vals, result_vals, result_weights);
    mkl_free(hamiltonian);
    mkl_free(sum_vals);
    mkl_free(set_flag);
    return nnz;
}

double perform_sampling(MKL_Complex16 *state, qaoa_data_t *meta_spec) {
    MKL_INT *choice_values;
    double result, best;
    double *sum_probs;
    double *cumul_probs;
    double *probabilities;
    double *temp;
    double *choices;
    MKL_INT space_dimension = meta_spec->machine_spec->space_dimension;

    probabilities = mkl_malloc(space_dimension * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(probabilities);
    sum_probs = mkl_malloc(meta_spec->cost_data->cx_range * sizeof(double), DEF_ALIGNMENT);
    check_alloc(sum_probs);
    cumul_probs = mkl_malloc(meta_spec->cost_data->cx_range * sizeof(double), DEF_ALIGNMENT);
    check_alloc(cumul_probs);
    choice_values = mkl_calloc((size_t) meta_spec->cost_data->cx_range, sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(choice_values);

    temp = mkl_calloc((size_t) meta_spec->run_spec->num_samples, sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(temp);
    choices = mkl_calloc((size_t) meta_spec->run_spec->num_samples, sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(choices);

    for (int i = 0; i < meta_spec->run_spec->num_samples; ++i) {
        temp[i] = 1.0;
    }
    compute_probabilities(state, meta_spec, probabilities);
    MKL_INT nnz = build_interval_array(meta_spec, probabilities, choice_values, sum_probs);

    cumul_probs[0] = sum_probs[0];
    for (int i = 1; i < nnz; ++i) {
        cumul_probs[i] = cumul_probs[i - 1] + sum_probs[i];
    }

    weighted_choice(choice_values, cumul_probs, nnz, meta_spec->run_spec, choices);
    best = choices[cblas_idamax(meta_spec->run_spec->num_samples, choices, 1)];
    if (best > meta_spec->qaoa_statistics->best_result) {
        meta_spec->qaoa_statistics->best_result = best;
    }

    mkl_free(probabilities);
    mkl_free(sum_probs);
    mkl_free(cumul_probs);
    mkl_free(choice_values);
    mkl_free(choices);
    mkl_free(temp);
    return best;
}

double evolve(unsigned num_params, const double *x, double *grad, qaoa_data_t *meta_spec){
    double result;
    int P = meta_spec->machine_spec->P;
    //Generate new initial state
    MKL_Complex16 *state = mkl_calloc((size_t)meta_spec->machine_spec->space_dimension, sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(state);
    initialise_state(state, meta_spec->machine_spec);
    //Apply our QAOA generation
    for(int i = 0; i < num_params/2; ++i){
        spmatrix_expm_z_diag(meta_spec->uc, x[i], meta_spec->machine_spec->space_dimension, state);
        spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16){x[i+P], 0.0}, (MKL_Complex16){0.0, -meta_spec->machine_spec->num_qubits},(MKL_Complex16){0.0, meta_spec->machine_spec->num_qubits}, meta_spec->machine_spec->space_dimension);
    }
    meta_spec->qaoa_statistics->num_evals++;
    //measure
    if(meta_spec->run_spec->sampling){
        result = perform_sampling(state, meta_spec);
    } else {
        result = expectation_value(state, meta_spec);
    }
    //teardown
    mkl_free(state);
    //Return single value;
    if(result > meta_spec->qaoa_statistics->best_result){
        meta_spec->qaoa_statistics->best_result = result;
    }
    if (meta_spec->run_spec->verbose) {
        printf("%d \n", meta_spec->qaoa_statistics->num_evals);
    }
    return result;
}

double evolve_restricted(unsigned num_params, const double *x, double *grad, qaoa_data_t *meta_spec) {
    double result;
    int P = meta_spec->machine_spec->P;
    //Generate new initial state
    MKL_Complex16 *state = mkl_calloc((size_t) meta_spec->machine_spec->space_dimension, sizeof(MKL_Complex16),
                                      DEF_ALIGNMENT);
    check_alloc(state);
    initialise_state(state, meta_spec->machine_spec);
    //Apply our QAOA generation
    for (int i = 0; i < (num_params - 1) / 2; ++i) {
        spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16) {x[i + P], 0.0},
                            (MKL_Complex16) {0.0, -meta_spec->machine_spec->num_qubits},
                            (MKL_Complex16) {0.0, meta_spec->machine_spec->num_qubits},
                            meta_spec->machine_spec->space_dimension);
        spmatrix_expm_z_diag(meta_spec->uc, x[i], meta_spec->machine_spec->space_dimension, state);
    }
    spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16) {x[num_params], 0.0},
                        (MKL_Complex16) {0.0, -meta_spec->machine_spec->num_qubits},
                        (MKL_Complex16) {0.0, meta_spec->machine_spec->num_qubits},
                        meta_spec->machine_spec->space_dimension);

    meta_spec->qaoa_statistics->num_evals++;
    //measure
    if (meta_spec->run_spec->sampling) {
        result = perform_sampling(state, meta_spec);
    } else {
        result = expectation_value(state, meta_spec);
    }
    //teardown
    mkl_free(state);
    //Return single value;
    if (result > meta_spec->qaoa_statistics->best_result) {
        meta_spec->qaoa_statistics->best_result = result;
    }
    if (meta_spec->run_spec->verbose) {
        printf("%d \n", meta_spec->qaoa_statistics->num_evals);
    }
    return result;
}
