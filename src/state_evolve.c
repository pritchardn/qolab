#include <mkl.h>
#include <mathimf.h>
#include "state_evolve.h"
#include "matrix_expm.h"

/**
 * @brief Returns the original problem hamiltonian as double values
 * @param uc The representation of the problem hamiltonian
 * @param hamiltonian The buffer to hold the final values
 * @param space_dimension The number of candidate solutions considered
 */
void extract_hamiltonian_double(MKL_Complex16 *uc, double *hamiltonian, MKL_INT space_dimension){
    for (int i = 0; i < space_dimension; ++i) {
        hamiltonian[i] = -uc[i].imag;
    }
}

/**
 * @brief Initializes a state vector as an equal superposition of all bit-strings
 * @param state The vector which will be initialised
 * @param mach_spec Data structure containing the number of qubits requested
 */
void initialise_state(MKL_Complex16 *state, machine_spec_t *mach_spec) {
    MKL_Complex16 init_value;
    init_value.real = 1.0 / sqrt(mach_spec->space_dimension);
    for(MKL_INT i = 0; i < mach_spec->space_dimension; ++i){
        state[i].real = init_value.real;
    }
}

/**
 * @brief Takes a complex state-vector and computes the complex conjugate; determining measurement probabilities
 * @param state The complex state taken as input
 * @param meta_spec The data structure containing simulation information
 * @param output The double array which will hold the result
 */
void compute_probabilities(MKL_Complex16 * state, qaoa_data_t *meta_spec, double *output){
    MKL_Complex16 *z_probabilities = mkl_malloc(meta_spec->machine_spec->space_dimension * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(z_probabilities);
    vzMulByConj(meta_spec->machine_spec->space_dimension, state, state, z_probabilities);
    vzAbs(meta_spec->machine_spec->space_dimension, z_probabilities, output);
    mkl_free(z_probabilities);
}

/**
 * @brief Checks that a given state vector is normalised
 * @param state The input state
 * @param meta_spec Used to determine the size of the state vector
 */
void check_probabilities(MKL_Complex16 *state, qaoa_data_t *meta_spec) {
    //TODO: Unit test for normalisation
    double result = 0.0;
    cblas_zdotc_sub(meta_spec->machine_spec->space_dimension, state, 1, state, 1, &result);
    printf("%f\n", result);
}

/**
 * @brief Determines the expectation value of measurement with respect to the problem Hamiltonian
 * @param state The complex state-vector
 * @param meta_spec The data-structure containing relevant information
 * @return Double value which is the expectation value of measurement
 */
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

/**
 * @brief Performs a binary search on a provided array for a particular target
 * @param array The target array to be searched
 * @param target The value to be searched for
 * @param size The size of the target array
 * @return The index of the hopefully found value
 */
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

/**
 * @brief Performs a weighted random choice from vals with respect to the provided weights
 * @param vals The output values to select from
 * @param weights Their corresponding probabilities
 * @param size The number of values to select from
 * @param run_spec Data structure containing the number of samples to perform
 * @param result The value selected
 * @warning Assumes that the values in vals and weights are correlated
 */
void weighted_choice(const MKL_INT *vals, const double *weights, MKL_INT size, run_spec_t *run_spec, double *result) {
    MKL_INT temp;
    double target;
    for (int i = 0; i < run_spec->num_samples; ++i) {
        target = rand() / (double) RAND_MAX;
        temp = binary_search(weights, target, size);
        result[i] = (vals[temp]);
    }
}

/**
 * @brief Builds the weights array required for a weighed choice
 * @param meta_spec Data structure containing simulation information
 * @param probabilities The array which holds the measurement probabilities for a given state
 * @param result_vals Will hold a packed vector containing the possible values of measurement
 * @param result_weights Will hold a packaed vector containing the corresponding sum probabilities of
 * measuring said values
 * @return The number of discrete result values detected
 */
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

/**
 * @brief Samples a provided state vector a number of times, recording the best value detected
 * @param state A complex state-vector
 * @param meta_spec Data structure containing simulation information and will hold resulting values
 * @return The best value detected while sampling
 */
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

    size_t index = cblas_idamax(meta_spec->cost_data->cx_range, sum_probs, 1);
    best = choice_values[index];
    if (best >= meta_spec->qaoa_statistics->best_result) {
        meta_spec->qaoa_statistics->best_result = best;
        if (sum_probs[index] > meta_spec->qaoa_statistics->best_result_prob) {
            meta_spec->qaoa_statistics->best_result_prob = sum_probs[index];
        }
    }

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

/**
 * @brief Performs a standard QAOA iteration (UBUC...)
 * @details Conforms to nlopt standards
 * @param num_params The number of optimimzation parameters present (2*P)
 * @param x The current candidate parameters
 * @param grad The gradient of the optimisation landscape (assumed to be NULL)
 * @param meta_spec Data structure containing all simulation information
 * @return Either the expectation value or sampled output value
 * @warning Assumes order of gamma(UC) then beta(UB) parameters.
 */
double evolve(unsigned num_params, const double *x, double *grad, qaoa_data_t *meta_spec){
    double result;
    int P = meta_spec->machine_spec->P;
    //Generate new initial state
    MKL_Complex16 *state = mkl_calloc((size_t)meta_spec->machine_spec->space_dimension, sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(state);
    initialise_state(state, meta_spec->machine_spec);
    //Apply our QAOA generation
    for(int i = 0; i < num_params / 2; ++i){
        spmatrix_expm_z_diag(meta_spec->uc, x[i], meta_spec->machine_spec->space_dimension, state);
        spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16){x[i+P], 0.0}, (MKL_Complex16){0.0, -meta_spec->machine_spec->num_qubits},(MKL_Complex16){0.0, meta_spec->machine_spec->num_qubits}, meta_spec->machine_spec->space_dimension);
    }
    meta_spec->qaoa_statistics->num_evals++;
    //measure
    if(meta_spec->run_spec->sampling){
        result = perform_sampling(state, meta_spec);
    } else {
        result = expectation_value(state, meta_spec);
        perform_sampling(state, meta_spec);
    }
    //teardown
    mkl_free(state);
    //Return single value;
    if (result > meta_spec->qaoa_statistics->best_expectation) {
        meta_spec->qaoa_statistics->best_expectation = result;
    }
    if (meta_spec->run_spec->verbose) {
        printf("%d \n", meta_spec->qaoa_statistics->num_evals);
    }
    return result;
}

/**
 * @brief Peroform a restricted QAOA iteration (UBUCUB...)
 * @details Simlar to a standard QAOA iteration but reverses the application of operators and applies one extra
 * UB operation. Conforms to nlopt standards
 * @param num_params The number of optimization parameters present (2*P + 1)
 * @param x The current candidate parameters
 * @param grad The gradient of the optimisation landscape (assumed to be NULL)
 * @param meta_spec Data structure containing all simulation information
 * @return Either the expectation value or sampled output value
 * @warning Assumes order of gamma(UC) then beta(UB) parameters.
 */
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
        perform_sampling(state, meta_spec);
    }
    //teardown
    mkl_free(state);
    //Return single value;
    if (result > meta_spec->qaoa_statistics->best_expectation) {
        meta_spec->qaoa_statistics->best_expectation = result;
    }
    if (meta_spec->run_spec->verbose) {
        printf("%d \n", meta_spec->qaoa_statistics->num_evals);
    }
    return result;
}
