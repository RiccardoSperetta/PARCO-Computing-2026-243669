#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

typedef uint32_t node_t;
typedef uint64_t edge_t;

typedef struct {
    node_t n_vertices;      // Number of vertices in graph
    edge_t n_edges;         // Number of edges in graph
    edge_t *row_ptr;        // Offsets into col_idx array (size: n_vertices + 1)
    node_t *col_idx;        // Destination vertices (size: n_edges)
} CSR;

typedef struct {
    node_t src;  // Source vertex
    node_t dst;  // Destination vertex
} Edge;

extern MPI_Datatype MPI_EDGE;
void init_mpi_datatypes(void);

//bfs utility:
int validate(CSR *g, node_t source, node_t local_start, node_t global_nverts, int* distance, int* parent);

//general utility:
void free_csr(CSR *csr);
void print_csr_stats(CSR *csr);


//parallel quick sort functions
int edge_compare(const void *a, const void *b);
long partition(Edge *arr, long low, long high);
void parallel_quicksort(Edge *arr, long low, long high, int depth);
void sort_edges(Edge *edges, long n_edges);



#endif
