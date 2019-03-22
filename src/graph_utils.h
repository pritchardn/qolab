/**
 * @author Nicholas Pritchard
 * @date 28/06/2018
 */

#ifndef GRAPHSIMILARITY_GRAPHUTIL_H
#define GRAPHSIMILARITY_GRAPHUTIL_H

#include <malloc.h>
#include <mkl.h>
#include "globals.h"

void print_graph(cost_data_t *cost_data, FILE *out);

void generate_graph(MKL_INT *graph1, int graph_size, float prob);

void copy_graph(const MKL_INT *src, MKL_INT *dest, int graph_size);

void generate_random(MKL_INT *graph, int graph_size, float prob);
void generate_undirected(MKL_INT *graph1, int graph_size, float prob);

void deform_flipped(const MKL_INT *src, MKL_INT *dest, int graph_size);
void deform_add_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int n);
void deform_add_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int n);
void deform_rem_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int n);
void deform_rem_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int n);

#endif //GRAPHSIMILARITY_GRAPHUTILS_H