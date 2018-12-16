//
// Created by nicholas on 14/12/18.
//

#include "reporting.h"

void machine_report(machine_spec_t * mach_spec, FILE *outfile){
    fprintf(outfile, "Machine Specification:\n"
                     "Qubits: %d\n"
                     "Decomposition: %d\n"
                     "Space Dimension: %lld\n", mach_spec->num_qubits, mach_spec->P, mach_spec->space_dimension);
}

void timing_report(qaoa_statistics_t *statistics, FILE *outfile){
    fprintf(outfile, "Timing report(s):\n"
                     "%f Total\n"
                     "%f UC\n"
                     "%f UB\n"
                     "%f Optimisation",
                     statistics->endTimes[0] - statistics->startTimes[0],
                     statistics->endTimes[1] - statistics->startTimes[1],
                     statistics->endTimes[2] - statistics->startTimes[2],
                     statistics->endTimes[3] - statistics->startTimes[3]);
}

void result_report(qaoa_statistics_t *statistics, FILE *outfile){

}

void optimiser_report(optimisation_spec_t *opt_spec, FILE *outfile){

}

void final_report(qaoa_data_t *meta_spec){
    machine_report(meta_spec->machine_spec, stdout);
    timing_report(meta_spec->qaoa_statistics, stdout);
}
