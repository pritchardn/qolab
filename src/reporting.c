//
// Created by nicholas on 14/12/18.
//

#include "reporting.h"

void machine_report(machine_spec_t * mach_spec, FILE *outfile){
    fprintf(outfile, "Machine Specification:\n"
                     "Qubits: %d\n"
                     "Decomposition: %d\n"
                     "Space Dimension: %lld", mach_spec->num_qubits, mach_spec->P, mach_spec->space_dimension);
}

void timing_report(qaoa_statistics_t *statistics, FILE *outfile){

}

void result_report(qaoa_statistics_t *statistics, FILE *outfile){

}

void optimiser_report(optimisation_spec_t *opt_spec, FILE *outfile){

}

void final_report(qaoa_data_t *meta_spec){
    printf("here\n");
}
