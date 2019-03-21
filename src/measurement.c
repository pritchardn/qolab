/**
 * @author Nicholas Pritchard
 * @date 20/03/2019
 */

#include "measurement.h"

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
 * @brief Returns the original problem hamiltonian as double values
 * @param uc The representation of the problem hamiltonian
 * @param hamiltonian The buffer to hold the final values
 * @param space_dimension The number of candidate solutions considered
 */
void extract_hamiltonian_double(MKL_Complex16 *uc, double *hamiltonian, MKL_INT space_dimension) {
    for (int i = 0; i < space_dimension; ++i) {
        hamiltonian[i] = -uc[i].imag;
    }
}

MKL_INT compact_probabilities(const double *probabilities, const double *hamiltonian, MKL_INT *values,
                              double *prob_compact, qaoa_data_t *meta_data) {
    MKL_INT nnz = 0;
    double *sum_vals;
    bool *set_flag;

    size_t num_vals = (size_t) meta_data->cost_data->cx_range;

    sum_vals = mkl_calloc(num_vals, sizeof(double), DEF_ALIGNMENT);
    set_flag = mkl_calloc(num_vals, sizeof(bool), DEF_ALIGNMENT);
    check_alloc(sum_vals);
    check_alloc(set_flag);

    for (MKL_INT i = 0; i < meta_data->machine_spec->space_dimension; ++i) {
        if (!set_flag[(int) hamiltonian[i] - 1]) {
            values[nnz] = (MKL_INT) hamiltonian[i];
            set_flag[(MKL_INT) hamiltonian[i] - 1] = true;
            nnz++;
        }
        sum_vals[(MKL_INT) hamiltonian[i] - 1] += probabilities[i];
    }

    vdPackV(nnz, sum_vals, values, prob_compact);

    mkl_free(sum_vals);
    mkl_free(set_flag);
    return nnz;
}

void cumulate_probabilities(const double *probabilities, double *cumul_probs, MKL_INT nnz) {
    cumul_probs[0] = probabilities[0];
    for (int i = 0; i < nnz; ++i) {
        cumul_probs[i] = cumul_probs[i - 1] + probabilities[i];
    }
}

MKL_INT weighted_choices(const MKL_INT *vals, const double *weights, MKL_INT nnz, run_spec_t *run_spec,
                         MKL_INT *result) {
    MKL_INT temp;
    MKL_INT sum = 0;
    double target;
    for (int i = 0; i < run_spec->num_samples; ++i) {
        target = rand() / (double) RAND_MAX;
        temp = binary_search(weights, target, nnz);
        result[i] = vals[temp];
        sum += result[i];
    }
    return sum;
}

/**
 * @brief Samples the provided probability distribution returning the best value sampled
 * @param probabilities
 * @param meta_spec
 * @return
 */
double sample(double *probabilities, qaoa_data_t *meta_spec) {
    double expectation;
    double *hamiltonian = NULL;
    double *prob_compact = NULL;
    double *cumul_probs = NULL;
    MKL_INT nnz;
    MKL_INT sample_sum;
    MKL_INT best_sample;
    MKL_INT *vals = NULL;
    MKL_INT *samples = NULL;


    MKL_INT space_dimension = meta_spec->machine_spec->space_dimension;
    vals = mkl_malloc(meta_spec->cost_data->cx_range * sizeof(MKL_INT), DEF_ALIGNMENT);
    prob_compact = mkl_calloc((size_t) meta_spec->cost_data->cx_range, sizeof(double), DEF_ALIGNMENT);

    check_alloc(vals);
    check_alloc(prob_compact);

    hamiltonian = mkl_malloc(space_dimension * sizeof(double), DEF_ALIGNMENT);
    extract_hamiltonian_double(meta_spec->uc, hamiltonian, space_dimension);

    nnz = compact_probabilities(probabilities, hamiltonian, vals, prob_compact, meta_spec);
    mkl_free(hamiltonian);
    cumul_probs = mkl_calloc((size_t) nnz, sizeof(double), DEF_ALIGNMENT);

    cumulate_probabilities(prob_compact, cumul_probs, nnz);
    mkl_free(prob_compact);

    //Allocate samples
    samples = mkl_malloc(meta_spec->run_spec->num_samples * sizeof(MKL_INT), DEF_ALIGNMENT);
    //Perform weighted samples
    sample_sum = weighted_choices(vals, cumul_probs, nnz, meta_spec->run_spec, samples);
    expectation = (double) sample_sum / (double) meta_spec->run_spec->num_samples;
    //Extract useful data
    best_sample = samples[cblas_idamax(meta_spec->run_spec->num_samples, (double *) samples, 1)];

    if (best_sample > meta_spec->qaoa_statistics->best_sample) {
        meta_spec->qaoa_statistics->best_sample = best_sample;
    }

    mkl_free(vals);
    mkl_free(cumul_probs);
    return expectation;
}

/**
 * @brief Determines the expectation value of measurement with respect to the problem Hamiltonian
 * @param state The complex state-vector
 * @param meta_data The data-structure containing relevant information
 * @return Double value which is the expectation value of measurement
 */
double expectation_value(double *probabilities, qaoa_data_t *meta_data) {
    MKL_INT space_dimension = meta_data->machine_spec->space_dimension;
    double *hamiltonian = NULL;
    double expectation;

    hamiltonian = mkl_malloc(space_dimension * sizeof(double), DEF_ALIGNMENT);
    extract_hamiltonian_double(meta_data->uc, hamiltonian, space_dimension);

    expectation = cblas_ddot(space_dimension, probabilities, 1, hamiltonian, 1);
    mkl_free(hamiltonian);

    return expectation;
}