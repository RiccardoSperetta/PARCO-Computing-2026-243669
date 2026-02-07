#include <mpi.h>
#include <omp.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include "graph_utils.h"
#include "bitset.h"

#define INF INT_MAX

char output_path[512];

int block_size, rank, p;
double start_time, local_sol_time;
double comm_time_start, local_comm_time = 0.0;

/* ============================================================================================================================
============================================================== BFS ============================================================
============================================================================================================================ */
int distributed_bfs(
    CSR g,                      // CSR graph
    node_t local_start,         // first global vertex owned
    node_t n,                   // number of total vertices
    node_t source,              // where the bfs has to start
    MPI_Comm comm
) {
    // Needed for normalizing row_ptr start-finish intervals
    node_t local_edges = g.row_ptr[0];

    edge_t traversed_edges = 0;
    //number of errors during BFS execution: all fine if = 0
    int local_result = 0;

    // Visited nodes tracked with a bitset: includes all nodes in the graph but is updated locally
    // = simplest compromise for ghost vertices management
    bitset_t visited = bitset_create(n);

    // Local distance vector - stores shortest path distance for each vertex from source
    atomic_int *d = aligned_alloc(64, g.n_vertices * sizeof(atomic_int));
    for (int i = 0; i < g.n_vertices; i++)
        d[i] = INF;

    // Local parent array - stores parent of each vertex in BFS tree
    int *parent = malloc(g.n_vertices * sizeof(int));
    for (int i = 0; i < g.n_vertices; i++)
        parent[i] = -1;

    // Local frontier 
    bitset_t frontier = bitset_create(g.n_vertices);
    node_t frontier_size = 0;

    local_comm_time = 0.0;
    start_time = MPI_Wtime();
    // Initialize: rank owning the source sets distance to 0 and adds to frontier
    if (source >= local_start && source < local_start + g.n_vertices) {
        node_t u = source - local_start;
        bitset_set(&visited, source);
        d[u] = 0;
        parent[u] = source;              // Source is its own parent
        bitset_set(&frontier, u);
        frontier_size = 1;
    }

    int level = 1;

    // Send buffers: each process prepares data for all other processes
    int *send_counts = calloc(p, sizeof(int));
    Edge **send_buf = malloc(p * sizeof(Edge*));


    /* ==============================================================
    main BFS loop
    ============================================================== */
    while (1) {
        edge_t tmp_edges = 0; //support variable for computing traversed edges

        /* ==============================================================
        * Exploration of local frontier with allocation for the outgoing edges of frontier nodes
        * ============================================================ */
    //beside from a trivial reduction and a small manual reduction on the send_counts array the parallelization is trivial
        #pragma omp parallel reduction(+:tmp_edges)
        {
            //if (rank == 0) printf("using %d threads\n", omp_get_num_threads());
            int *local_send = calloc(p, sizeof(int));

            #pragma omp for schedule(guided, 128) nowait
            for (node_t i = 0; i < g.n_vertices; i++) {
                if (!bitset_test(&frontier, i)) continue;

                tmp_edges += g.row_ptr[i+1] - g.row_ptr[i];

                for (edge_t e = g.row_ptr[i] - local_edges; e < g.row_ptr[i+1] - local_edges; e++) {
                    node_t v = g.col_idx[e];
                    if (bitset_test(&visited, v)) continue;
                    int owner = v / block_size;
                    local_send[owner]++;
                }
            }

            #pragma omp critical
            for (int i = 0; i < p; i++) {
                send_counts[i] += local_send[i];
            }

            free(local_send);
        }
        
        // Allocate send_buf based on calculated send_counts -> avoiding over allocation of n_edges for all ranks...
        // still an over-estimation: we may count multiple times the same vertex
        for (int i = 0; i < p; i++) {
            send_buf[i] = malloc(send_counts[i] * sizeof(Edge));
            send_counts[i] = 0; // Reset send_counts for reuse in the next loop 
        }

        // Populate send_buf with edges to be sent
    // Parallelizing this loop appeared to be not worthy because of the multiple contentions on
    // visited bitset but most importantly the allocation of send_buf which would have to be subject of a reduction (one for each rank)
        for (node_t i = 0; i < g.n_vertices; i++) {
            if(bitset_test(&frontier, i)) {
                for (edge_t e = g.row_ptr[i] - local_edges; e < g.row_ptr[i + 1] - local_edges; e++) {
                    node_t v = g.col_idx[e];            // neighbor of node in frontier (global)
                    if (bitset_test(&visited, v) == 1) {continue;}
                    bitset_set(&visited, v);

                    int owner = v / block_size;         // process owning the destination vertex of the edge
                    
                    send_buf[owner][send_counts[owner]].src = i+local_start;
                    send_buf[owner][send_counts[owner]++].dst = v;
                }
            }
        }
        tmp_edges *= 2; // it's just the same scan as before considering this time I'm storing the edges I'm going to send
        traversed_edges += tmp_edges;

        /* ==============================================================
        * Exchange of NUMBER of outgoing edges before Alltoallv -> each rank knows how much edges it needs to expect
        * ============================================================ */
        int *recv_counts = malloc(p * sizeof(int));
        comm_time_start = MPI_Wtime();
        MPI_Alltoall(send_counts, 1, MPI_INT,
                     recv_counts, 1, MPI_INT,
                     comm);
        local_comm_time += (MPI_Wtime() - comm_time_start);
        // => each rank knows how many "external visits" it's going to receive


        // Calculate displacements for send/receive buffers
        // needed for a correct distribution (and receivement) in the Alltoallv between all the p processes
        int *sdispls = malloc(p * sizeof(int));     //SEND
        int *rdispls = malloc(p * sizeof(int));     //RECEIVE

        sdispls[0] = rdispls[0] = 0;
        for (int i = 1; i < p; i++) {
            sdispls[i] = sdispls[i - 1] + send_counts[i - 1];
            rdispls[i] = rdispls[i - 1] + recv_counts[i - 1];
        }

        int total_recv = rdispls[p - 1] + recv_counts[p - 1];    // = total "external visists" received
        Edge *recv_buf = malloc(total_recv * sizeof(Edge));

        // Flatten send buffer from 2D to 1D
        int total_send = sdispls[p - 1] + send_counts[p - 1];    // = total "external visits" that must be comunicated to the original owner
        Edge *send_flat = malloc(total_send * sizeof(Edge));

        for (int i = 0; i < p; i++)
            for (int j = 0; j < send_counts[i]; j++)
                send_flat[sdispls[i] + j] = send_buf[i][j];

        
        /* ==============================================================
        * Exchange of ACTUAL outgoing edges before Alltoallv -> all ranks can now update their frontier
        * ============================================================ */
        comm_time_start = MPI_Wtime();
        MPI_Alltoallv(send_flat, send_counts, sdispls, MPI_EDGE,
                      recv_buf, recv_counts, rdispls, MPI_EDGE,
                      comm);
        local_comm_time += (MPI_Wtime() - comm_time_start);

    
        // Build next frontier: add unvisited vertices and set parent pointers
        frontier_size = 0;
        bitset_clear(&frontier);
        
        int th_frontier_size = 0;
    // the parellalization here is trivial if not for the race conditions on both distance vector and the frontier
    // which can be atomic, implying high contention on high degree vertices
        #pragma omp parallel for reduction(+:th_frontier_size) schedule(static)
        for (int i = 0; i < total_recv; i++) {
            Edge e = recv_buf[i];
            if (e.dst < local_start || e.dst >= local_start + g.n_vertices) {
                local_result++; //RACE CONDITION - not a real problem I just need it to be >0 if something went wrong
                continue;
            }
            node_t local_v = e.dst - local_start;
            int expected = INF;

            //instead of going for an omp critical section:
            if (atomic_compare_exchange_strong(&d[local_v], &expected, level)) {
                //only ONE thread will access this code
                parent[local_v] = e.src;
                bitset_set_atomic(&frontier, local_v); // = avoids race conditions on the update of bitset frontier
                th_frontier_size++; 
            }
        }

        frontier_size += th_frontier_size;

        /* ==============================================================
        * Check global termination: stop if no new vertices were discovered = global frontier is empty
        * ============================================================== */ 
        node_t global_frontier;
        comm_time_start = MPI_Wtime();
        MPI_Allreduce(&frontier_size, &global_frontier, 1,
                      MPI_UINT32_T, MPI_SUM, comm);
        local_comm_time += (MPI_Wtime() - comm_time_start);

        if (global_frontier == 0) {
            local_sol_time = MPI_Wtime() - start_time; //rank individual solution time (they should all be quite close)
        }

        // Cleanup temporary buffers
        for (int i = 0; i < p; i++)
            free(send_buf[i]);
        free(recv_counts);
        free(sdispls);
        free(rdispls);
        free(send_flat);
        free(recv_buf);

        if (global_frontier == 0) { // = all processes don't have any more vertices to explore
#ifdef DEBUG
            if(rank == 0) printf("BFS finished at level %d\n", level);   
            printf("RANK %d - traversed %ld edges\n", rank, traversed_edges);
#endif
            break;
        }

        level++;

    }

    /* ==============================================================
    * Post BFS validation
    * ============================================================== */
    local_result += validate(&g, source, local_start, n, (int*)d, parent);

#ifdef DEBUG
    if (local_result == 0) {
        printf("All fine on rank %d\n", rank);
    } else {
        printf("Something's wrong in rank %d: %d errors found\n", rank, local_result);
    }
#endif
        
    int total_result = 0;
    MPI_Reduce(&local_result, &total_result, 1, MPI_INT, MPI_SUM, 0, comm);


    /* ==============================================================
    * Printing measurements on the result file
    * ============================================================== */

    double total_time, total_comm_time;
    double TEPS;
    double max_over_mean, cv = 0;
#ifdef DEBUG
    printf("RANK %d - traversed %lu edges\n", rank, traversed_edges);
#endif
    compute_imbalance_metrics(local_sol_time, local_comm_time, traversed_edges, &total_time, &total_comm_time, &TEPS, &max_over_mean, &cv, comm);

    if (rank == 0) { //rank responsible for writing results
        FILE* results = fopen(output_path, "a");
        if (results == NULL) {
            perror("Error opening file");
            return 1;
        }
        fprintf(results, "%e %e %e %f %f\n", total_time, total_comm_time, TEPS, max_over_mean, cv);
    }

    free(send_buf);
    free(send_counts);
    bitset_free(&visited);
    bitset_free(&frontier);
    free(d);
    free(parent);

    return total_result;
}

/* ============================================================================================================================
============================================================== MAIN ===========================================================
============================================================================================================================ */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    init_mpi_datatypes();

    if (argc != 3) {
        printf("Usage: %s <csr_file> <run_specs>\n", argv[0]);
        return 1;
    }

    if (rank == 0) printf("using %d OpenMP threads per rank\n", omp_get_max_threads());

    const char *filename = argv[1];
    const char *run_specs = argv[2];
    const char *last_slash = strrchr(filename, '/');
    const char *basename = last_slash ? last_slash + 1 : filename;
    char *dot = strrchr(basename, '.');
    int basename_len = dot ? (dot - basename) : strlen(basename);

    if (strncmp(basename, "kronecker", 9) == 0) {
        snprintf(output_path, sizeof(output_path), "results/weak_scaling/%s.txt", run_specs);
    } else {
        snprintf(output_path, sizeof(output_path), "results/%.*s/%s.txt", basename_len, basename, run_specs);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    CSR graph;
    graph.row_ptr = NULL;
    graph.col_idx = NULL;
    graph.n_vertices = 0;
    graph.n_edges = 0;

    node_t total_vertices;

    /* ============================================================== 
    STEP 1: each rank reads its portion of CSR through MPI IO
    ============================================================== */
    MPI_File handle;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &handle);

    // manually making the node_t and edge_t types match with MPI datatypes
    MPI_File_read_at(handle, 0, &total_vertices, 1, MPI_UINT32_T, NULL);

    srand(11);
    node_t search_keys[64];
    for (int i = 0; i < 64; i++) {
        search_keys[i] = rand() % total_vertices; //source is always randomized across all possible graph IDs
    }

    block_size = (total_vertices + p - 1) / p; //ceil division
    // now each rank can compute its offset in order to start reading row_ptr:

    long long base_offset = sizeof(node_t) + sizeof(edge_t);
    long long row_offset = base_offset + (block_size * rank) * sizeof(edge_t);
    graph.n_vertices = block_size;      //each process will own a "block_size" of vertices
    if (block_size * (rank + 1) > total_vertices) {                   //we want to avoid reading too many vertices in the last process, which may get fewer nodes
        graph.n_vertices = total_vertices - block_size * (rank);
    }

    graph.row_ptr = malloc((graph.n_vertices+1) * sizeof(edge_t));
    
    MPI_File_read_at(handle, row_offset, graph.row_ptr,   
                    graph.n_vertices+1, MPI_UINT64_T, NULL);
    // now each rank also knows what portion of col_idx it must take care of:

    edge_t start, end;
    start = graph.row_ptr[0];
    end = graph.row_ptr[graph.n_vertices];
    graph.n_edges = end-start;          //each process will own all the outgoing edges of the nodes they own
    graph.col_idx = malloc(graph.n_edges * sizeof(node_t));

    long long col_offset = sizeof(node_t) + sizeof(edge_t) + (total_vertices+1) * sizeof(edge_t) + (start) * sizeof(node_t);

    MPI_File_read_at(handle, col_offset, graph.col_idx,   
                graph.n_edges, MPI_UINT32_T, NULL);

    MPI_File_close(&handle);

    #ifdef DEBUG
    printf("RANK %d - %d vertices and %ld edges, going from %ld to %ld\n", rank, graph.n_vertices, graph.n_edges, start, start + graph.n_edges);
    print_csr_stats(&graph);
    #endif

    /* ============================================================== 
    STEP 2: BFS algorithm
    ============================================================== */
    int result = 0;
    for (int i=0; i<64; i++) {
#ifdef DEBUG
        if (rank == 0) printf("run %d\n", i);
#endif
        result += distributed_bfs(graph, rank*block_size, total_vertices, search_keys[i], MPI_COMM_WORLD);
        if (result != 0) {
            fprintf(stderr, "BFS #%d produced an incorrect BFS tree\n", i);
            break;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    

    free_csr(&graph);

    MPI_Type_free(&MPI_EDGE);
    MPI_Finalize();
    if(result != 0) { //something, based on the validate function, went wrong
        fprintf(stderr, "Something went wrong with this BFS\n");
        return 1;
    }

    return 0;
}