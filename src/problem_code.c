#include <stdio.h>
#include "problem_code.h"

/**
 * To be implemented by the user. This defines the problem investigated by the QAOA
 * @param i
 * @param num_qubits
 * @param cost_data
 * @return A single integer value
 */

#define CHECK_BIT(var, pos) ((var) & (1 << (pos)))

int Cx(int i, int num_qubits, cost_data_t *cost_data){
    int result = 0;
    for (int j = 0; j < num_qubits; ++j) {
        if (CHECK_BIT(i, j)) {
            for (int k = 0; k < num_qubits; ++k) {
                if (!CHECK_BIT(i, k)) {
                    if (cost_data->graph[j * num_qubits + k] == 1) {
                        result++;
                    }
                }
            }
        } else {
            for (int k = 0; k < num_qubits; ++k) {
                if (CHECK_BIT(i, k)) {
                    if (cost_data->graph[j * num_qubits + k] == 1) {
                        result++;
                    }
                }
            }
        }
    }
    printf("%d %d\n", i, result);
    return result;
}

bool mask(unsigned int i, cost_data_t *cost_data) {
    return true;
}