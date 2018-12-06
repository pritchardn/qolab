#include <mkl.h>
#include <stdbool.h>

typedef struct {
    double startTimes[5];
    double endTimes[5];
    double result, classical_exp, random_exp;
    int max_value, max_index;
    int term_status;
} qaoa_statistics;

typedef struct {
    bool timing;
    bool correct;
    bool report;
} run_spec_t;

typedef struct {
    int num_qubits;
    int P;
    MKL_INT space_dimension;
} machine_spec_t;

typedef struct {
    int nlopt_method;
    double xtol;
    double ftol;
} optimisation_spec_t;


typedef struct{
    machine_spec_t *machine_spec;
    run_spec_t *run_spec;
}qaoa_data_t;

void qaoa(machine_spec_t *mach_spec, optimisation_spec_t *opt_spec, run_spec_t *run_spec);
