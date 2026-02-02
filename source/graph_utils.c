#include "graph_utils.h"
#include "mpi.h"
#include <omp.h>
#include <stdio.h>

// TODO: better commenting

MPI_Datatype MPI_EDGE;  // definition (storage) here

void init_mpi_datatypes(void) {
    MPI_Type_contiguous(2, MPI_UINT32_T, &MPI_EDGE);   // change MPI_INT if node_t is different
    MPI_Type_commit(&MPI_EDGE);
}

int validate(CSR *g, node_t source, node_t local_start, node_t global_nverts, int* distance, int* parent) {
    int local_errors = 0;
    edge_t local_edges = g->row_ptr[0];

    // TODO: maybe consider only some vertices
    for (node_t v = 0; v < g->n_vertices; v++) {
        node_t global_v = v + local_start;           
        node_t p = parent[v];

        if (p == -1) continue;                    // unreachable

        if (global_v == source) {                      // source must be parent of itself and be at distance = 0
            if (p != source || distance[global_v] != 0) local_errors++;
        } else {
            if (p < 0 || p >= global_nverts || p == global_v) { // out of range condition + multiple sources
                local_errors++;
                continue;
            }

            // does edge exist? - works ONLY for UNDIRECTED graphs!
            int edge_exists = 0;
            for (edge_t j = g->row_ptr[v]-local_edges; j < g->row_ptr[v+1]-local_edges; j++) {
                if (g->col_idx[j] == p) {
                    edge_exists = 1;
                    break;
                }
            }
            if (edge_exists != 1) {
                local_errors++;
                continue;
            }

            // BFS tree distance invariant
            if (p >= local_start && p < local_start + g->n_vertices) { //for local vertices immediate check
                if (distance[v] != distance[p-local_start]+1) {
                    local_errors++;
                }
            } 
            //parent may be remote → defer to query
            // or skip for now if you don't want cross-rank yet
            else {

            }
        }
    }

    return local_errors;
}

void free_csr(CSR *g) {
    if (g) {
        free(g->row_ptr);
        free(g->col_idx);
    }
}

// Print basic statistics about the CSR
void print_csr_stats(CSR *g) {
    printf("\n=== CSR Statistics ===\n");
    printf("Vertices: %d\n", g->n_vertices);
    printf("Edges: %ld\n", g->n_edges);
    printf("Average degree: %.2f\n", (double)g->n_edges / g->n_vertices);
    
    // Find min/max degree
    node_t min_degree = g->n_edges;
    node_t max_degree = 0;
    for (int i = 0; i < g->n_vertices; i++) {
        long degree = g->row_ptr[i + 1] - g->row_ptr[i];
        if (degree < min_degree) min_degree = degree;
        if (degree > max_degree) max_degree = degree;
    }
    printf("Min degree: %d\n", min_degree);
    printf("Max degree: %d\n", max_degree);

    printf("row_ptr: [");
    node_t cap = 15 < g->n_vertices ? 15: g->n_vertices;
    for (int i=0; i<=cap; i++) {
        printf("%ld ", g->row_ptr[i]);
    }
    printf("]\n");
    printf("col_idx: [");
    node_t cap2 = 30 < g->n_edges ? 30: g->n_edges;
    for (int i=0; i<cap2; i++) {
        printf("%d ", g->col_idx[i]);
    }
    printf("]\n");
}


/* ============================================================================
 * PARALLEL QUICKSORT IMPLEMENTATION
 * ============================================================================
 * For 2^20+ edges, serial qsort becomes a bottleneck.
 * This uses OpenMP tasks to parallelize the recursive quicksort calls.
 */

// Standard comparison for edges: sort by source, then by destination
int edge_compare(const void *a, const void *b) {
    Edge *ea = (Edge *)a;
    Edge *eb = (Edge *)b;
    if (ea->src != eb->src)
        return ea->src - eb->src;
    return ea->dst - eb->dst;
}

// Partition function for quicksort
long partition(Edge *arr, long low, long high) {
    Edge pivot = arr[high];
    long i = low - 1;
    
    for (long j = low; j < high; j++) {
        if (edge_compare(&arr[j], &pivot) <= 0) {
            i++;
            Edge temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    Edge temp = arr[i + 1];
    arr[i + 1] = arr[high];
    arr[high] = temp;
    return i + 1;
}

// Parallel quicksort using OpenMP tasks
void parallel_quicksort(Edge *arr, long low, long high, int depth) {
    if (low < high) {
        long pi = partition(arr, low, high);
        
        // Only create parallel tasks for larger chunks and reasonable depth
        // This avoids task overhead for small subarrays
        if (depth > 0 && (high - low) > 10000) {
            #pragma omp task shared(arr)
            parallel_quicksort(arr, low, pi - 1, depth - 1);
            
            #pragma omp task shared(arr)
            parallel_quicksort(arr, pi + 1, high, depth - 1);
            
            #pragma omp taskwait
        } else {
            // Fall back to serial sort for small chunks
            parallel_quicksort(arr, low, pi - 1, 0);
            parallel_quicksort(arr, pi + 1, high, 0);
        }
    }
}

// Wrapper to start parallel sort
void sort_edges(Edge *edges, long n_edges) {
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Depth limits parallelism to avoid too many tasks
            // log2(num_threads) * 2 is a reasonable heuristic
            int max_depth = omp_get_num_threads() > 1 ? 
                           __builtin_ctz(omp_get_num_threads()) * 2 : 0;
            parallel_quicksort(edges, 0, n_edges - 1, max_depth);
        }
    }
}
