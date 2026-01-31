#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "graph_utils.h"

// TODO: better commenting
CSR* edgelist_to_csr(const char *filename, int directed) {
    printf("Converting edge list to CSR format...\n");
    double start_time = omp_get_wtime();
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }
    
    /* ========================================================================
     * STEP 1: Count edges and find maximum vertex ID
     * ======================================================================== 
     * since SNAP graphs may not be contiguous we need to first read the graph 
     * and find out what's the max ID -> then we'll be able to map those IDs from 0 to n-1
     * = we avoid dealing with extrimely sparse graphs simply because of the labels used for their vertices
     */
    printf("STEP 1: Scanning file...\n");
    node_t max_vertex = 0;
    edge_t edge_count = 0;
    node_t u, v;
    char line[256];
    
    while (fgets(line, sizeof(line), fp)) {
        // Skip comment lines (common in SNAP datasets)
        if (line[0] == '#') continue;
        
        // Try to parse two integers (source and destination)
        if (sscanf(line, "%d %d", &u, &v) == 2) {
            edge_count++;
            // Track the maximum vertex ID we've seen
            if (u > max_vertex) max_vertex = u;
            if (v > max_vertex) max_vertex = v;
        }
    }
        
    // For undirected graphs, we store each edge twice (u->v and v->u)
    edge_t n_edges = directed ? edge_count : 2 * edge_count;
    
    printf("Found %ld edges with vertices with a maximum ID of %d\n", n_edges, max_vertex);
    
    // We're sure on the number of edges of the graph:
    Edge *edges = malloc(n_edges * sizeof(Edge));
    if (!edges) {
        fprintf(stderr, "Error: Cannot allocate memory for edges\n");
        fclose(fp);
        return NULL;
    }
    
    /* ========================================================================
     * STEP 2: Read all edges into memory + mapping of indices to [0...n-1]
     * ======================================================================== 
     * We rewind the file and read again, this time storing edges.
     */
    printf("STEP 2: storing edge list...\n");
    rewind(fp);
    edge_t idx = 0;


    /* allocate map */
    node_t *map = malloc((max_vertex + 1) * sizeof(node_t));
    for (node_t i = 0; i <= max_vertex; i++) {
        map[i] = -1;
    }
    
    /* remap edges as soon as you copy them in the edges array */
    node_t next_id = 0;

    //TODO: first shuffle, THEN relabel in order to allow for basic load balancing!
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        
        if (sscanf(line, "%d %d", &u, &v) == 2) {
            //update mapping if new vertex is found
            if (map[u] == -1) map[u] = next_id++;
            u = map[u];
            if (map[v] == -1) map[v] = next_id++;
            v = map[v];

            // Store the edge
            edges[idx].src = u;
            edges[idx].dst = v;
            idx++;
            
            // For undirected graphs, also store the reverse edge
            if (!directed) {
                edges[idx].src = v;
                edges[idx].dst = u;
                idx++;
            }
        }
    }
    fclose(fp);
    
    /* ========================================================================
     * STEP 3: SORT EDGES by source vertex (and destination as tiebreaker)
     * ======================================================================== 
     * CSR format requires edges from the same source to be 
     * contiguous in memory. After sorting:
     */
    printf("STEP 3: Sorting %ld edges...\n", n_edges);
#ifdef _OPENMP
    double sort_start = omp_get_wtime();
    sort_edges(edges, n_edges);
    //qsort(edges, n_edges, sizeof(Edge), edge_compare);
    printf("Sorting took %.2f seconds\n", omp_get_wtime() - sort_start);
#endif

#ifdef DEBUG
    int cap = 20 > n_edges ? n_edges : 20;
    for (int i=0; i<cap; i++) {
        printf("%d %d\n", edges[i].src, edges[i].dst);
    }
#endif

    /* ========================================================================
     * STEP 4: BUILD CSR STRUCTURE
     * ======================================================================== 
     * CSR has two arrays:
     * 
     * 1. row_ptr[i] = index in col_idx where vertex i's edges start
     *    - Size: n_vertices + 1
     *    - row_ptr[n_vertices] = total number of edges
     *    - Edges for vertex i are in col_idx[row_ptr[i] ... row_ptr[i+1]-1]
     * 
     * 2. col_idx[j] = destination vertex for edge j
     *    - Size: n_edges
     *    - Stores all destination vertices in sorted order
     */
    
    //TODO: also save this edge list in data/edge_list/filename.txt (?)
    // needed for graph500
    printf("Building CSR arrays...\n");
    CSR *csr = malloc(sizeof(CSR));
    if (!csr) {
        free(edges);
        return NULL;
    }
    
    printf("TOTAL vertices found = %d\n", next_id);
    node_t n_vertices = next_id;    //only by mapping all unique vertices we can understand how many nodes the graph really has
    csr->n_vertices = n_vertices; 
    csr->n_edges = n_edges;
    
    // Allocate row_ptr (initialized to 0)
    csr->row_ptr = calloc(n_vertices + 1, sizeof(edge_t));
    csr->col_idx = malloc(n_edges * sizeof(node_t));
    
    if (!csr->row_ptr || !csr->col_idx) {
        fprintf(stderr, "Error: Cannot allocate CSR arrays\n");
        free(edges);
        free_csr(csr);
        return NULL;
    }
    
    // TODO: it's easy to do this in one step
    /* ========================================================================
     * STEP 1: Fill col_idx and count edges per vertex
     * ======================================================================== 
     * Since edges are sorted by source, we can:
     * - Copy all destinations to col_idx
     * - Count how many edges each vertex has
     * 
     * We use row_ptr[src+1] temporarily to count, will fix in next step.
     */
    for (edge_t i = 0; i < n_edges; i++) {
        csr->col_idx[i] = edges[i].dst;
        // Count edges for this source vertex (store in next position)
        csr->row_ptr[edges[i].src + 1]++;
    }
    
    /* ========================================================================
     * STEP 2: Convert counts to offsets via prefix sum
     * ======================================================================== 
     * Transform row_ptr from "counts" to "cumulative offsets"
     * 
     * Before: row_ptr = [0, 2, 1, 2] (counts for vertices 0,1,2)
     * After:  row_ptr = [0, 2, 3, 5] (starting positions)
     * 
     * This is a prefix sum (cumulative sum) operation.
     */
    for (node_t i = 1; i <= n_vertices; i++) {
        csr->row_ptr[i] += csr->row_ptr[i - 1];
    }
    
    // Free the temporary edge array - we don't need it anymore
    free(edges);
    
    printf("Memory usage: row_ptr=%.2f MB, col_idx=%.2f MB\n",
           (n_vertices + 1) * sizeof(edge_t) / 1024.0 / 1024.0,
           n_edges * sizeof(node_t) / 1024.0 / 1024.0);
    
    return csr;
}


/* ========================================================================================
MAIN
======================================================================================== */

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Usage: %s <edge_list_file>\n", argv[0]);
        return 1;
    }
    
    // Set number of OpenMP threads (or use OMP_NUM_THREADS env variable)
    // omp_set_num_threads(16);
    
    printf("Using %d OpenMP threads\n", omp_get_max_threads());
    
    // Convert edge list to CSR (1 = directed, 0 = undirected)
    CSR *graph = edgelist_to_csr(argv[1], 1);
    
    if (graph) {
        #ifdef DEBUG
        print_csr_stats(graph);
        #endif

        // Extract filename from path and create output path
        char output_path[512];
        const char *filename = argv[1];
        const char *last_slash = strrchr(filename, '/'); // = returns a pointer to the last occurrence of the character c in the string s.
        const char *basename = last_slash ? last_slash + 1 : filename;
        char *dot = strchr(basename, '.');
        int basename_len = dot ? (dot - basename) : strlen(basename);

        snprintf(output_path, sizeof(output_path), "data/csr/%.*s.bin", basename_len, basename);

        // Write CSR to binary file
        FILE *out_fp = fopen(output_path, "wb");
        if (!out_fp) {
            fprintf(stderr, "Error: Cannot open output file %s\n", output_path);
        } else {
            fwrite(&graph->n_vertices, sizeof(node_t), 1, out_fp);
            fwrite(&graph->n_edges, sizeof(edge_t), 1, out_fp);
            fwrite(graph->row_ptr, sizeof(edge_t), graph->n_vertices + 1, out_fp);
            fwrite(graph->col_idx, sizeof(node_t), graph->n_edges, out_fp);
            fclose(out_fp);
            printf("CSR saved to %s\n", output_path);
        }

        free_csr(graph);
    }
    
    return 0;
}

/* 
what's going on:
    run program with data/raw/file_name.txt 
    first read of file in order to understand what's the max-vertex ID and number of edges
    load edges in memory and meanwhile re-label vertices
    sort (even parallely if needed) 
    build CSR
    print (if wanted) 
    eventually store on data/csr/file_name.bin
*/