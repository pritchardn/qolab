#ifndef QOLAB_STATE_EVOLVE_H
#define QOLAB_STATE_EVOLVE_H
#include "globals.h"

double evolve(unsigned num_params, const double *x, double *grad, qaoa_data_t *meta_spec);

#endif //QOLAB_STATE_EVOLVE_H
