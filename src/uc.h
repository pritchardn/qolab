//
// Created by nicholas on 28/06/18.
//

#ifndef GRAPHSIMILARITY_UC_H
#define GRAPHSIMILARITY_UC_H

#include <stdio.h>
#include <mathimf.h>
#include <complex.h>
#include <mkl.h>
#include <limits.h>
#include "globals.h"
#include "qaoa.h"
#include "problem_code.h"

void generate_uc(int num_qubits, MKL_Complex16 *uc, MKL_INT (*Cx)(MKL_INT, int, cost_data *), cost_data *info,
                 qaoa_statistics *stats);
#endif //GRAPHSIMILARITY_UC_H
