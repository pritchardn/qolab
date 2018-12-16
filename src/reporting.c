//
// Created by nicholas on 14/12/18.
//

#include <time.h>
#include "reporting.h"

FILE *file_generate(qaoa_data_t *meta_spec){
    FILE *rtn;
    struct tm *now;
    time_t seconds;
    char buffer[128];

    seconds = time(NULL);
    now = localtime(&seconds);
    buffer[0] = '\0';
    sprintf(buffer, "Q%dP%dM%dT%d%d%d%d%d.out",
            meta_spec->machine_spec->num_qubits,
            meta_spec->machine_spec->P,
            meta_spec->opt_spec->nlopt_method,
            now->tm_mday, now->tm_mon+1, now->tm_hour, now->tm_min, now->tm_sec);

    rtn = fopen(buffer, "w+");

    if(rtn==NULL){
        perror("Attempting to open report file");
    }
    printf("Writing to%s\n",buffer);
    return rtn;
}

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
                     "%f Optimisation\n"
                     "%d #Evals\n",
                     statistics->endTimes[0] - statistics->startTimes[0],
                     statistics->endTimes[1] - statistics->startTimes[1],
                     statistics->endTimes[2] - statistics->startTimes[2],
                     statistics->endTimes[3] - statistics->startTimes[3],
                     statistics->num_evals);
}

void result_report(qaoa_statistics_t *statistics, FILE *outfile){
    fprintf(outfile, "Result report:\n"
                     "%d %d gOpt, Loc\n"
                     "%f Final Result\n"
                     "%f Best Result\n"
                     "%f Classical Exp\n"
                     "%f Initial Exp\n"
                     "%d Termination value\n",
                     statistics->max_value, statistics->max_index,
                     statistics->result,
                     statistics->best_result,
                     statistics->classical_exp,
                     statistics->random_exp,
                     statistics->term_status);
}

void optimiser_report(optimisation_spec_t *opt_spec, int P, FILE *outfile){
    fprintf(outfile, "Optimisation report\n"
                     "%d Method\n"
                     "%d Max evals\n"
                     "%f xtol\n"
                     "%f ftol\n",
                     opt_spec->nlopt_method,
                     opt_spec->max_evals,
                     opt_spec->xtol,
                     opt_spec->ftol);
    int i,j;
    for (i = 0; i < P; ++i) {
        fprintf(outfile, "%f ", opt_spec->parameters[i]);
    }fprintf(outfile, "Gammas\n");
    for (j = 0; j < P; ++j) {
        fprintf(outfile, "%f ", opt_spec->parameters[j+P]);
    }fprintf(outfile, "Betas\n");
}

void final_report(qaoa_data_t *meta_spec){
    FILE *oFile;
    if(meta_spec->run_spec->report){
        oFile = file_generate(meta_spec);
        meta_spec->run_spec->outfile = oFile;
    } else {
        meta_spec->run_spec->outfile = stdout;
    }
    machine_report(meta_spec->machine_spec, meta_spec->run_spec->outfile);
    timing_report(meta_spec->qaoa_statistics, meta_spec->run_spec->outfile);
    optimiser_report(meta_spec->opt_spec, meta_spec->machine_spec->P, meta_spec->run_spec->outfile);
    result_report(meta_spec->qaoa_statistics, meta_spec->run_spec->outfile);
    if(meta_spec->run_spec->report){
        fprintf(meta_spec->run_spec->outfile, "Test\n");
        fclose(meta_spec->run_spec->outfile);
    }
}
