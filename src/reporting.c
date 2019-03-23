/**
 * @author Nicholas Pritchard
 * @date 1/07/2018
 * @brief Contains methods used to print out reports
 */
#include <time.h>
#include "reporting.h"

/**
 * @brief Generates a filename for a given run
 * @param meta_spec Contains specifications about the simulation
 * @return A new file pointer
 */
FILE *file_generate(qaoa_data_t *meta_spec){
    time_t seconds;
    FILE *result;
    struct tm *now;
    char buffer[128];

    seconds = time(NULL);
    now = localtime(&seconds);
    buffer[0] = '\0';
    sprintf(buffer, "Q%dP%dM%dT%d%d%d%d%d.out",
            meta_spec->machine_spec->num_qubits,
            meta_spec->machine_spec->P,
            meta_spec->opt_spec->nlopt_method,
            now->tm_mday, now->tm_mon + 1, now->tm_hour, now->tm_min, now->tm_sec);

    result = fopen(buffer, "w+");

    if (result == NULL) {
        perror("Attempting to open report file");
    }
    printf("Writing to%s\n", buffer);
    return result;
}

/**
 * @brief Reports on the machine specification
 * @param mach_spec Contains the machine specification
 * @param outfile The file stream to print to
 */
void machine_report(machine_spec_t *mach_spec, FILE *outfile) {
    if (outfile == NULL) {
        outfile = stdout;
    }
    fprintf(outfile, "Machine Specification:\n"
                     "Qubits: %d\n"
                     "Decomposition: %d\n"
                     "Space Dimension: %d\n", mach_spec->num_qubits, mach_spec->P, mach_spec->space_dimension);
}

/**
 * @brief Reports on the timing of the simulation
 * @param statistics Contains the timing information
 * @param outfile The file stream to print to
 */
void timing_report(qaoa_statistics_t *statistics, FILE *outfile){
    if (outfile == NULL) {
        outfile = stdout;
    }
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

/**
 * @brief Reports on the correctness of the algorithm
 * @param statistics Contains the correctness information
 * @param outfile The file stream to print to
 */
void result_report(qaoa_statistics_t *statistics, FILE *outfile){
    if (outfile == NULL) {
        outfile = stdout;
    }
    fprintf(outfile, "Result report:\n"
                     "%d %d gOpt, Loc\n"
                     "%f Final Expectation\n"
                     "%d Best Sample\n"
                     "%f Best expectation\n"
                     "%f Classical Exp\n"
                     "%f Initial Exp\n",
            statistics->max_value, statistics->max_index,
            statistics->result,
            (MKL_INT) statistics->best_sample,
            statistics->best_expectation,
            statistics->classical_exp,
            statistics->random_exp);
    nlopt_termination_parser(statistics->term_status, outfile);
}

/**
 * @brief Reports on the final state of the classical optimisation
 * @param opt_spec Contains the informaiton about the classical optimiser
 * @param P The amount of decomposition used in the simulation
 * @param outfile The file stream to print to
 */
void optimiser_report(optimization_spec_t *opt_spec, int P, FILE *outfile) {
    if (outfile == NULL) {
        outfile = stdout;
    }
    fprintf(outfile, "Optimisation report\n"
                     "%d Method\n"
                     "%d Max evals\n"
                     "%.2e xtol\n"
                     "%.2e ftol\n",
            opt_spec->nlopt_method,
            opt_spec->max_evals,
            opt_spec->xtol,
            opt_spec->ftol);
    int i,j;
    for (i = 0; i < P; ++i) {
        fprintf(outfile, "%f ", opt_spec->parameters[i]);
    }fprintf(outfile, "Gammas\n");
    for (j = 0; j < P; ++j) {
        fprintf(outfile, "%f ", opt_spec->parameters[j + P]);
    }fprintf(outfile, "Betas\n");
}

/**
 * @brief Reports on an individual optimisation iteration
 * @param measurement The most recent measurement value
 * @param meta_spec Contains all information about the simulation
 */
void iteration_report(double measurement, qaoa_data_t *meta_spec) {
    FILE *oFile = meta_spec->run_spec->outfile;
    if (meta_spec->run_spec->outfile == NULL) {
        oFile = stdout;
    }
    fprintf(oFile, "%d: %f\n", meta_spec->qaoa_statistics->num_evals, measurement);
}

/**
 * @brief Print a full report after the simulation has been run
 * @param meta_spec Contains information about the simulation
 */
void final_report(qaoa_data_t *meta_spec){
    FILE *oFile;
    if(meta_spec->run_spec->report){
        oFile = file_generate(meta_spec);
        meta_spec->run_spec->outfile = oFile;
    } else {
        meta_spec->run_spec->outfile = stdout;
    }
    machine_report(meta_spec->machine_spec, meta_spec->run_spec->outfile);
    if (meta_spec->run_spec->timing) {
        timing_report(meta_spec->qaoa_statistics, meta_spec->run_spec->outfile);
    }
    if (meta_spec->run_spec->correct) {
        optimiser_report(meta_spec->opt_spec, meta_spec->machine_spec->P, meta_spec->run_spec->outfile);
        result_report(meta_spec->qaoa_statistics, meta_spec->run_spec->outfile);
    }
    if(meta_spec->run_spec->report){
        fclose(meta_spec->run_spec->outfile);
    }
}

/**
 * @brief Parses the nlopt termination code for a sensible printout
 * @param nlopt_code The code to be parsed
 * @param outfile The file stream to print to
 */
void nlopt_termination_parser(nlopt_result nlopt_code, FILE *outfile) {
    if (outfile == NULL) {
        outfile = stderr;
    }
    switch (nlopt_code) {
        case NLOPT_SUCCESS:
            fprintf(outfile, "1 nlopt_success\n");
            break;
        case NLOPT_STOPVAL_REACHED:
            fprintf(outfile, "2 nlopt_stopval_reached\n");
            break;
        case NLOPT_FTOL_REACHED:
            fprintf(outfile, "3 nlopt_ftol_reached\n");
            break;
        case NLOPT_XTOL_REACHED:
            fprintf(outfile, "4 nlopt_xtol_reached\n");
            break;
        case NLOPT_MAXEVAL_REACHED:
            fprintf(outfile, "5 nlopt_maxeval_reached\n");
            break;
        case NLOPT_MAXTIME_REACHED:
            fprintf(outfile, "6 nlopt_maxtime_reached\n");
            break;
        case NLOPT_FAILURE:
            fprintf(outfile, "-1 nlopt_failure\n");
            break;
        case NLOPT_INVALID_ARGS:
            fprintf(outfile, "-2 nlopt_invalid_args\n");
            break;
        case NLOPT_OUT_OF_MEMORY:
            fprintf(outfile, "-3 nlopt_out_of_memory\n");
            break;
        case NLOPT_ROUNDOFF_LIMITED:
            fprintf(outfile, "-4 nlopt_roundoff_limited\n");
            break;
        case NLOPT_FORCED_STOP:
            fprintf(outfile, "-5 nlopt_forced_stop\n");
            break;
        default:
            fprintf(outfile, "%d nlopt_other\n", nlopt_code);
            break;
    }
}
