//
// Created by nicholas on 28/06/18.
//
#include "uc.h"

//TODO: Use machine_spec structs
void generate_uc(int num_qubits, MKL_Complex16 *uc, MKL_INT (*Cx)(MKL_INT, int, cost_data *), cost_data *info,
                 qaoa_statistics *stats) {
    int current;
    MKL_INT i;
    MKL_INT space_dimension = (MKL_INT) pow(2, num_qubits);
    double c_sum = 0.0, classic_prob;
    classic_prob = (double) 1.0 / factorial(info->x_range);
    stats->max_value = INT_MIN;
    stats->max_index = 0;
    for (i = 0; i < space_dimension; ++i) {
        current = Cx(i, num_qubits, info);
        if (current > stats->max_value) {
            stats->max_value = current;
            stats->max_index = i;
        }
        c_sum += (double) current * classic_prob;
        //Pre-multplying by -I to save time later.
        uc[i].imag = -current;
    }
    stats->classical_exp = c_sum;
    stats->random_exp = c_sum * 1 / classic_prob / pow(2, num_qubits);
}
