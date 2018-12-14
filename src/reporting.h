//
// Created by nicholas on 14/12/18.
//

#ifndef QOLAB_REPORTING_H
#define QOLAB_REPORTING_H
#include "globals.h"

void machine_report(machine_spec_t * mach_spec);
void timing_report(qaoa_statistics_t *statistics);
void result_report(qaoa_statistics_t *statistics);
void optimiser_report(optimisation_spec_t *opt_spec);

void iteration_report(qaoa_data_t *meta_spec);

void final_report(qaoa_data_t *meta_spec);

#endif //QOLAB_REPORTING_H
