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

void generate_uc(qaoa_data_t *meta_data, cost_data_t *cost_data, int (*Cx)(int, int, cost_data_t *));

#endif //GRAPHSIMILARITY_UC_H
