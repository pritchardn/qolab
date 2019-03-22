/**
 * @author Nicholas Pritchard
 * @date 28/06/28
 * @brief A number of utility functions for generating and manipulating graphs
 */

#include <time.h>
#include "graph_utils.h"

/**
 * @brief Generates an array of random doubles
 * @param num_request The number of random doubles to be generated
 * @param buffer The buffer (assumed to be at least size num_request) to be filled with numbers
 * @warning Uses srand, should move to a new generator soon
 */
void random_doubles(int num_request, double *buffer) {
    srand((unsigned) time(0)); //TODO: Change random engine to something reasonable
    for (int i = 0; i < num_request; ++i) {
        buffer[i] = rand() / (double) RAND_MAX;
    }
}

/**
 * @brief Pretty-print for our graph data-structure
 * @details Assumes a 1D array representing a flattened adjacency matrix.
 * @param cost_data The meta-structure which should contain a graph
 * @param out The filestream to print to
 */
void print_graph(cost_data_t *cost_data, FILE *out){
    int graph_size = cost_data->num_vertices;
    for (int i = 0; i < graph_size; ++i) {
        for (int j = 0; j < graph_size; ++j) {
            fprintf(out, "%d ", cost_data->graph[graph_size * i + j]);
        }
        fprintf(out, "\n");
    }fprintf(out, "\n");
}

/**
 * @brief Generates two random graphs. Generates two ErdosRenyi graphs using system random. Should be avoided due to well
 * known limitations of system rand.
 * @param graph1  Assumed to be of graph_size * graph_size in length
 * @param graph2  Assumed to be of graph_size * graph_size in length
 * @param graph_size The number of vertices required
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


/**
 * @brief Generates a random erdos-renyi graph
 * @details Generates a random number for each edge, populates those which are < prob
 * @param graph The adjacency matrix to hold the graph
 * @param graph_size The number of vertices in the graph
 * @param prob The population probability of each edge
 */
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


/**
 * @brief Generates an undirected graph by edge population
 * @details Generates random numbers for the upper triangle of the graph adjacency matrix. Those < prob are populated
 * along with their mirrored entry (a[i][j] and a[j][i] are both set).
 * @param graph1 The graph to be populated
 * @param graph_size The number of vertices in the graph
 * @param prob The probability of each edge begin populated
 */
void generate_undirected(MKL_INT *graph1, int graph_size, float prob) {
    double *randomNums = mkl_calloc((size_t) factorial(graph_size), sizeof(double), DEF_ALIGNMENT);
    check_alloc(randomNums);
    random_doubles(factorial(graph_size), randomNums);
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

/**
 * @brief Copies the adjacency matrix from one graph to another
 * @param src The original graph
 * @param dest The destination graph
 * @param graph_size The number of vertices in the src graph
 * @warning Assumes the size of dest is at least that of src
 */
void copy_graph(const MKL_INT *src, MKL_INT *dest, int graph_size) {
    for (int i = 0; i < graph_size * graph_size; ++i) {
        dest[i] = src[i];
    }
}


/**
 * @brief Counts the number of zero entries in a graph adjacency matrix
 * @param src The graph to be traversed
 * @param graph_size The number of vertices in the graph
 * @return The number of zero entries
 */
int count_0(const MKL_INT *src, int graph_size) {
    int count_0 = 0;
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (src[i] == 0) {
            count_0++;
        }
    }
    return count_0;
}

/**
 * @brief Counts the number of ones in a graph adjacency matrix
 * @param src The graph to be traversed
 * @param graph_size The number of vertices in the graph
 * @return The number of ones in the adjacency matrix
 */
int count_1(const MKL_INT *src, int graph_size) {
    int count_1 = 0;
    for (int i = 0; i < graph_size * graph_size; ++i) {
        if (src[i] == 1) {
            count_1++;
        }
    }
    return count_1;
}


/**
 * @brief Adds a random edge to a directed graph
 * @param src The original graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param count The number of edges to be added
 */
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

/**
 * @brief Adds a random edge to an undirected graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param count The number of edges to add
 */
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

/**
 * @brief Removes random edges from a directed graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param count The number of edges to remove
 */
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

/**
 * @brief Removes random edges from an undirected graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param count The number of edges to remove
 */
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


/**
 * @brief Flips the adjacency matrix of a graph in the vertical axis
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 */
void deform_flipped(const MKL_INT *src, MKL_INT *dest, int graph_size) {
    for (int i = 0; i < graph_size; ++i) {
        for (int j = 0; j < graph_size; ++j) {
            //Flip each entry along the vertical axis
            dest[graph_size * i + j] = src[graph_size * i + (graph_size + j)];
        }
    }
}


/**
 * @brief Adds directed edges to a graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param n The number of edges to add
 */
void deform_add_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_0(src, graph_size);
        add_edge_direct(src, dest, graph_size, count);
    }
}

/**
 * @brief Adds undirected edges to a graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param n The number of edges to add
 */
void deform_add_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_0(src, graph_size) / 2;
        add_edge_undirect(src, dest, graph_size, count);
    }
}

/**
 * @brief Removes directed edges from a graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param n The number of edges to remove
 */
void deform_rem_direct(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_1(src, graph_size);
        rem_edge_direct(src, dest, graph_size, count);
    }
}


/**
 * @brief Removes edges from an undirected graph
 * @param src The graph to be examined
 * @param dest The graph to be edited
 * @param graph_size The number of vertices in both graphs
 * @param n The number of edges to remove
 */
void deform_rem_undirect(const MKL_INT *src, MKL_INT *dest, int graph_size, int n) {
    int count;
    for (int i = 0; i < n; ++i) {
        count = count_1(src, graph_size) / 2;
        rem_edge_undirect(src, dest, graph_size, count);
    }
}


