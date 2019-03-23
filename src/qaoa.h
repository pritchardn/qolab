/**
 * @author Nicholas Pritchard
 * @date 1/07/18
 */

#include <mkl.h>
#include <stdbool.h>
#include <time.h>
#include "problem_code.h"
#include "globals.h"

void qaoa(machine_spec_t *mach_spec, cost_data_t *cost_data, optimization_spec_t *opt_spec, run_spec_t *run_spec,
          bool retain);
/*! \mainpage Main Menu
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Depedencies
 *
 * First one must make sure the following dependencies are met.
 *   - Intel MKL
 *   - Nlopt
 *   - MPI (if running on a cluster)
 *
 * \subsection step2 Step 2: Defining your problem
 *
 * The next step is to define your problem. The functions in problem_code.c need to be implemented and the struct in
 * problem_code.h may need to be adapted to contain additional information.
 *
 * For more elaboration refer to the appropriate documentation pages.
 *
 * \subsection step3 Step 3: Compiling
 *
 * We recommend using the intel compiler where possible however gcc will suffice.
 *
 * We recommend using a simple makefile over complex build-tools. This gets users more familiar with the process and is
 * in general simpler.
 *
 * \subsubsection desktop_build Desktop
 * In the \verbatim /build \endverbatim directory we provide an example makefile for a desktop build of Qolab.
 *
 * Please ensure your environment variables are set for the various libraries
 * \subsubsection hpc_build Cluster
 *
 * In the \verbatim /build_hpc \endverbatim director we provide an exampke makefile and slurm script for a clust build
 * of Qolab
 *
 * Again, ensure that the appropriate packages and environment variables are set.
 *
 * To deploy a successful job script, each compute node should be running a single MPI process with maximal threads
 * enabled.
 * \section wishlist_sec Wishlist
 *   - Windows support
 *
 * \section contact_sec Contact
 * Do not hesitate to reach out report errors and suggest features
 *   - Github: https://github.com/pritchardn/qolab
 *
 * \section licensing Licence
 * GPLv3
 */