#include "uc.h"

/**
 * \brief Generates the solution hamiltonian which encodes the problem dependent solutions to every possible bit-string.
 * @param meta_data Contains all the information about our simulation
 * @param Cx The function which implements the problem-dependent cost-function
 * @param mask (Optional) A bit-string mask (the same as UB-generation) to avoid computing the cost-function for
 * invalid candidate solutions in the restricted QAOA.
 * @warning Will probably hit double precision for the c_sum statistic very quickly
 */
void generate_uc(qaoa_data_t *meta_data, int (*Cx)(int, int, cost_data_t *),
                 bool (*mask)(unsigned int, cost_data_t *cost_data)) {
    int current, i, num_qubits;
    num_qubits = meta_data->machine_spec->num_qubits;
    double c_sum = 0.0, classic_prob;
    classic_prob = (double) 1.0 / meta_data->cost_data->x_range;
    meta_data->qaoa_statistics->max_value = INT_MIN;
    meta_data->qaoa_statistics->max_index = 0;
    for (i = 0; i < meta_data->machine_spec->space_dimension; ++i) {
        current = Cx(i, num_qubits, meta_data->cost_data);
        //current = 0;
        if (current > meta_data->qaoa_statistics->max_value && mask((unsigned) i, meta_data->cost_data)) {
            meta_data->qaoa_statistics->max_value = current;
            meta_data->qaoa_statistics->max_index = i;
        }
        c_sum += (double) current * classic_prob;
        //Pre-multplying by -I to save time later.
        meta_data->uc[i].imag = -current;
    }
    meta_data->qaoa_statistics->classical_exp = c_sum;
    meta_data->qaoa_statistics->random_exp = c_sum * 1 / classic_prob / pow(2, num_qubits);
}
