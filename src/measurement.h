/**
 * @author Nicholas Pritchard
 * @date 20/03/2019
 */

#ifndef QOLAB_SAMPLING_H
#define QOLAB_SAMPLING_H

#include <mkl.h>
#include "globals.h"

double sample(double *probabilities, qaoa_data_t *meta_spec);

double expectation_value(double *probabilities, qaoa_data_t *meta_spec);

#endif //QOLAB_SAMPLING_H
