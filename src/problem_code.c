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

/**
 * @brief An optional function to be implemented by the user. Determines whether a given candidate solution is valid
 * @param i The candidate solution
 * @param cost_data Contains problem dependent data, may or may not be useful
 * @return True if the given solution is valid, false otherwise
 */
bool mask(unsigned int i, cost_data_t *cost_data) {
    return true;
}