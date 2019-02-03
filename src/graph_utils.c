//
// Created by nicholas on 28/06/18.
//

#include <time.h>
#include "graph_utils.h"
#include "globals.h"

void random_doubles(int num_request, double *buffer) {
    srand((unsigned) time(0));
    for (int i = 0; i < num_request; ++i) {
        buffer[i] = rand() / (double) RAND_MAX;
    }
}

/**
 * Generates two random graphs. Generates two ErdosRenyi graphs using system random. Should be avoided due to well
 * known limitations of system rand.
 * @param graph1  Assumed to be of graph_size * graph_size in length
 * @param graph2  Assumed to be of graph_size * graph_size in length
 * @param graph_size
 * @param prob  Defaults to 0.5
 */
void generate_graph(MKL_INT *graph1, int graph_size, float prob) {
    double *randomNums = mkl_calloc((size_t) 2 * graph_size * graph_size, sizeof(double), DEF_ALIGNMENT);
    check_alloc(randomNums);
    random_doubles(graph_size * graph_size * 2, randomNums);
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (randomNums[i] < prob) {
            graph1[i] = 1;
        }
    }
    mkl_free(randomNums);
}


void generate_random(MKL_INT *graph, int graph_size, float prob) {
    double *randomNums = mkl_calloc((size_t) graph_size * graph_size, sizeof(double), DEF_ALIGNMENT);
    check_alloc(randomNums);
    random_doubles(graph_size * graph_size, randomNums);
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (randomNums[i] < prob) {
            graph[i] = 1;
        }
    }
    mkl_free(randomNums);
}


void generate_undirected(MKL_INT *graph1, int graph_size, float prob) {
    double *randomNums = mkl_calloc((size_t) factorial(graph_size), sizeof(double), DEF_ALIGNMENT);
    check_alloc(randomNums);
    random_doubles((int) factorial(graph_size), randomNums);
    for (int i = 0; i < graph_size; ++i) {
        for (int j = i + 1; j < graph_size; ++j) {
            if (randomNums[i * graph_size + j] < prob) {
                graph1[graph_size * i + j] = 1;
                graph1[graph_size * j + i] = 1;
            }
        }
    }
    mkl_free(randomNums);
}

void copy_graph(const MKL_INT *src, MKL_INT *dest, int graph_size) {
    for (int i = 0; i < graph_size * graph_size; ++i) {
        dest[i] = src[i];
    }
}


int count_0(const MKL_INT *src, MKL_INT *dest, int graph_size) {
    int count_0 = 0;
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (src[i] == 0) {
            count_0++;
        }
    }
    return count_0;
}


int count_1(const MKL_INT *src, MKL_INT *dest, int graph_size) {
    int count_1 = 0;
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (src[i] == 1) {
            count_1++;
        }
    }
    return count_1;
}


void add_edge_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int count) {
    int choice, count_z = 0;
    choice = rand() % count;
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (src[i] == 0) {
            if (choice == count_z) {
                dest[i] = 1;
                break;
            } else {
                count_z++;
            }
        }
    }
}

void add_edge_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int count) {
    int choice, count_z = 0;
    choice = rand() % count;
    for (int i = 0; i < graph_size; ++i) {
        for (int j = 0; j < graph_size; ++j) {
            if (src[i * graph_size + j] == 0) {
                if (choice == count_z) {
                    dest[i * graph_size + j] = 1;
                    dest[j * graph_size + i] = 1;
                    break;
                } else {
                    count_z++;
                }
            }
        }
    }
}


void rem_edge_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int count) {
    int choice, count_o = 0;
    choice = rand() % count;
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (src[i] == 1) {
            if (choice == count_o) {
                dest[i] = 0;
                break;
            } else {
                count_o++;
            }
        }
    }
}

void rem_edge_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int count) {
    int choice, count_0 = 0;
    choice = rand() % count;
    for (int i = 0; i < graph_size; ++i) {
        for (int j = 0; j < graph_size; ++j) {
            if (src[i * graph_size + j] == 1) {
                if (choice == count_0) {
                    dest[i * graph_size + j] = 0;
                    dest[j * graph_size + i] = 0;
                    break;
                } else {
                    count_0++;
                }
            }
        }
    }
}


void deform_flipped(const MKL_INT *src, MKL_INT *dest, int graph_size) {
    for (int i = 0; i < graph_size; ++i) {
        for (int j = 0; j < graph_size; ++j) {
            //Flip each entry along the vertical axis
            dest[graph_size * i + j] = src[graph_size * i + (graph_size + j)];
        }
    }
}


void deform_add_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_0(src, dest, graph_size);
        add_edge_direct(src, dest, graph_size, count);
    }
}

void deform_add_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_0(src, dest, graph_size) / 2;
        add_edge_undirect(src, dest, graph_size, count);
    }
}


void deform_rem_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_1(src, dest, graph_size);
        rem_edge_direct(src, dest, graph_size, count);
    }
}


void deform_rem_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_1(src, dest, graph_size) / 2;
        rem_edge_undirect(src, dest, graph_size, count);
    }
}


