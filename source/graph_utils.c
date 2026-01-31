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

void free_csr(CSR *csr) {
    if (csr) {
        free(csr->row_ptr);
        free(csr->col_idx);
    }
}

// Print basic statistics about the CSR
void print_csr_stats(CSR *csr) {
    printf("\n=== CSR Statistics ===\n");
    printf("Vertices: %d\n", csr->n_vertices);
    printf("Edges: %ld\n", csr->n_edges);
    printf("Average degree: %.2f\n", (double)csr->n_edges / csr->n_vertices);
    
    // Find min/max degree
    long min_degree = csr->n_edges;
    long max_degree = 0;
    for (int i = 0; i < csr->n_vertices; i++) {
        long degree = csr->row_ptr[i + 1] - csr->row_ptr[i];
        if (degree < min_degree) min_degree = degree;
        if (degree > max_degree) max_degree = degree;
    }
    printf("Min degree: %ld\n", min_degree);
    printf("Max degree: %ld\n", max_degree);

    printf("row_ptr: [");
    for (int i=0; i<=csr->n_vertices; i++) {
        printf("%ld ", csr->row_ptr[i]);
    }
    printf("]\n");
    printf("col_idx: [");
    for (int i=0; i<csr->n_edges; i++) {
        printf("%d ", csr->col_idx[i]);
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
