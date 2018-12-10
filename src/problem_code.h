//
// Created by nicholas on 6/12/18.
//
#ifndef QOLAB_PROBLEM_CODE_H
#define QOLAB_PROBLEM_CODE_H
#include "globals.h"

typedef struct {
    int *graph;
    MKL_INT x_range;
    MKL_INT cx_range;
} cost_data_t;

int Cx(int i, int num_qubits, cost_data_t *cost_data);

#endif //QOLAB_PROBLEM_CODE_H
