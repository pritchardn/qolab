#ifndef QOLAB_REPORTING_H
#define QOLAB_REPORTING_H
#include "globals.h"

void machine_report(machine_spec_t * mach_spec, FILE *outfile);
void timing_report(qaoa_statistics_t *statistics, FILE *outfile);
void result_report(qaoa_statistics_t *statistics, FILE *outfile);
void optimiser_report(optimisation_spec_t *opt_spec, int P, FILE *outfile);

void nlopt_termination_parser(nlopt_result nlopt_code, FILE *outfile);

void iteration_report(qaoa_data_t *meta_spec);

void final_report(qaoa_data_t *meta_spec);

#endif //QOLAB_REPORTING_H
