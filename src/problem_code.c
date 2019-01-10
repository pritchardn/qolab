#include "problem_code.h"

/**
 * To be implemented by the user. This defines the problem investigated by the QAOA
 * @param i
 * @param num_qubits
 * @param cost_data
 * @return A single integer value
 */
int Cx(int i, int num_qubits, cost_data_t *cost_data){
    return i;
}

bool mask(unsigned int i) {
    return true;
}