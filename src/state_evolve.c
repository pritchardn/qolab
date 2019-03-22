#include <mkl.h>
#include <mathimf.h>
#include "state_evolve.h"
#include "matrix_expm.h"
#include "measurement.h"

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
    if (result != 1.0) {
        printf("Not normalised\n");
        exit(EXIT_FAILURE);
    }
    //printf("%f\n", result);
}

double measure(MKL_Complex16 *state, qaoa_data_t *meta_data) {
    double result;
    double *probabilities = mkl_malloc(meta_data->machine_spec->space_dimension * sizeof(double), DEF_ALIGNMENT);

    check_alloc(probabilities);

    compute_probabilities(state, meta_data, probabilities);

    if (meta_data->run_spec->sampling) {
        //Perform sampling
        result = sample(probabilities, meta_data);
    } else {
        //Perform expectation value
        result = expectation_value(probabilities, meta_data);
    }

    mkl_free(probabilities);
    return result;
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
        spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16) {x[i + P], 0.0},
                            (MKL_Complex16) {0.0, -meta_spec->ub_eigenvalue},
                            (MKL_Complex16) {0.0, meta_spec->ub_eigenvalue}, meta_spec->machine_spec->space_dimension);
    }
    meta_spec->qaoa_statistics->num_evals++;
    //measure
    result = measure(state, meta_spec);
    //teardown
    mkl_free(state);
    //Return single value;
    if (result > meta_spec->qaoa_statistics->best_sample) {
        meta_spec->qaoa_statistics->best_sample = result;
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
    check_probabilities(state, meta_spec);
    //Apply our QAOA generation
    for (int i = 0; i < (num_params - 1) / 2; ++i) {
        spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16) {x[i + P], 0.0},
                            (MKL_Complex16) {0.0, -meta_spec->ub_eigenvalue},
                            (MKL_Complex16) {0.0, meta_spec->ub_eigenvalue},
                            meta_spec->machine_spec->space_dimension);
        spmatrix_expm_z_diag(meta_spec->uc, x[i], meta_spec->machine_spec->space_dimension, state);
    }
    spmatrix_expm_cheby(&meta_spec->ub, state, (MKL_Complex16) {x[num_params], 0.0},
                        (MKL_Complex16) {0.0, -meta_spec->ub_eigenvalue},
                        (MKL_Complex16) {0.0, meta_spec->ub_eigenvalue},
                        meta_spec->machine_spec->space_dimension);
    meta_spec->qaoa_statistics->num_evals++;
    //measure
    result = measure(state, meta_spec);
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
