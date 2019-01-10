//
// Created by nicholas on 1/07/18.
//

#ifndef GRAPHSIMILARITY_UB_H
#define GRAPHSIMILARITY_UB_H

#include <stdio.h>
#include <mathimf.h>
#include <complex.h>
#include <mkl.h>
#include <stdbool.h>
#include "globals.h"

void generate_ub(qaoa_data_t *meta_data, bool (*mask)(unsigned int));

#endif //GRAPHSIMILARITY_UB_H
