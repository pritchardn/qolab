#ifndef QOLAB_PROBLEM_CODE_H
#define QOLAB_PROBLEM_CODE_H

#include <mkl.h>
#include <stdbool.h>
/**
 * A struct which contains additional information for the cost function to work
 * x_range and cx_range must both be included for all cases.
 */
typedef struct {
    MKL_INT x_range;    // The upper bound of cost-function input values.
    MKL_INT cx_range;   // The upper bound of cost-function output values.
    MKL_INT *graph;
    MKL_INT num_vertices;
} cost_data_t;

int Cx(int i, int num_qubits, cost_data_t *cost_data);

bool mask(unsigned int i, cost_data_t *cost_data);

#endif //QOLAB_PROBLEM_CODE_H
