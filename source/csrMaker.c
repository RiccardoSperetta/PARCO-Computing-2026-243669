#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>
#include "graph_utils.h"

CSR* edgelist_to_csr(const char *filename, int directed, int shuffle) {      // REMEMBER: for later pull technique we MUST work with undirected graphs!
    printf("Converting edge list to CSR format...\n");
    
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
    
    // NOW we can be sure on the number of edges of the graph:
    Edge *edges = malloc(n_edges * sizeof(Edge));
    if (!edges) {
        fprintf(stderr, "Error: Cannot allocate memory for edges\n");
        fclose(fp);
        return NULL;
    }
    
    /* ========================================================================
     * STEP 2: Read all edges into memory 
     * ======================================================================== 
     */
    printf("STEP 2: storing edge list...\n");
    rewind(fp);
    edge_t idx = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        
        if (sscanf(line, "%d %d", &u, &v) == 2) {
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
     * STEP 3: SHUFFLING of edges for basic load balancing 
     * ======================================================================== 
     */
    
    printf("STEP 3: Shuffling edges...");
    if(shuffle) {
        srand(omp_get_wtime() + getpid());
        printf("\n");
        for (edge_t i = n_edges - 1; i > 0; i--) {
            edge_t j = rand() % (i + 1); // Pick a random index from 0 to i
            // Swap edges[i] and edges[j]
            Edge temp = edges[i];
            edges[i] = edges[j];
            edges[j] = temp;
        }
    } else {
        printf("skip\n");
    }
    


    /* ========================================================================
     * STEP 4: MAPPING of indices to [0...n-1]
     * ======================================================================== 
     * After shuffling, we impose a new remapping of those vertices so that we can 
     * first sort and then build a CSR accordingly to this new renaming of vertices
     */
    printf("STEP 4: Mapping vertices IDs...\n");
     /* allocate map */
    node_t *map = malloc((max_vertex + 1) * sizeof(node_t));
    for (node_t i = 0; i <= max_vertex; i++) {
        map[i] = -1;
    }
    node_t next_id = 0;

    for (edge_t i=0; i<n_edges; i++) {
        node_t u = edges[i].src;
        node_t v = edges[i].dst;
        //update mapping if new vertex is found
        if (map[u] == -1) map[u] = next_id++;
        edges[i].src = map[u];
        if (map[v] == -1) map[v] = next_id++;
        edges[i].dst = map[v];
    }

    
    /* ========================================================================
     * STEP 5: SORT EDGES by source vertex (and destination as tiebreaker)
     * ======================================================================== 
     * CSR format requires edges from the same source to be 
     * contiguous in memory, that's why we need orting:
     */
    printf("STEP 5: Sorting %ld edges...\n", n_edges);
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
     * STEP 6: BUILD CSR STRUCTURE
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
    
    printf("STEP 6: Building CSR arrays...\n");
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
    
    // most straight forward approach: first I compute the col_idx array
    for (edge_t i = 0; i < n_edges; i++) {
        csr->col_idx[i] = edges[i].dst;
        // Count edges for this source vertex (store in next position)
        csr->row_ptr[edges[i].src + 1]++;
    }

    //then I compute the prefix sum on row_ptr
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
    if (argc != 4) {
        printf("Usage: %s <edge_list_file> <directed> <shuffle>\n", argv[0]);
        return 1;
    }

    int directed = atoi(argv[2]);
    int shuffle = atoi(argv[3]);
    
    // Set number of OpenMP threads (or use OMP_NUM_THREADS env variable)
    // omp_set_num_threads(16);
    
    printf("Using %d OpenMP threads\n", omp_get_max_threads());
    
    // Convert edge list to CSR (1 = directed, 0 = undirected)
    CSR *graph = edgelist_to_csr(argv[1], directed, shuffle);
    
    if (graph) {
        #ifdef DEBUG
        print_csr_stats(graph);
        #endif

        // Extract filename from path and create output path
        char output_path[512];
        const char *filename = argv[1];
        const char *last_slash = strrchr(filename, '/'); // = returns a pointer to the last occurrence of the character c in the string s.
        const char *basename = last_slash ? last_slash + 1 : filename;
        char *dot = strrchr(basename, '.');
        int basename_len = dot ? (dot - basename) : strlen(basename);

        if (shuffle) {
            snprintf(output_path, sizeof(output_path), "data/csr/%.*s_shuffled.bin", basename_len, basename);
        } else {
            snprintf(output_path, sizeof(output_path), "data/csr/%.*s.bin", basename_len, basename);
        }

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
