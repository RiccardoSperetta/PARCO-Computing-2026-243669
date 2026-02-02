#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include "graph_utils.h"
#include "bitset.h"

#define INF UINT_MAX

int block_size, rank, p;
double start_time, end_time;
double comm_time_start, comm_time_end;

/* ============================================================================================================================
============================================================== BFS ============================================================
============================================================================================================================ */
void distributed_bfs(
    CSR g,                      // CSR graph
    node_t local_start,         // first global vertex owned
    node_t n,                   // number of total vertices
    node_t source,              // where the bfs has to start
    MPI_Comm comm
) {
    // Needed for normalizing row_ptr start-finish intervals
    node_t local_edges = g.row_ptr[0];
    // Visited nodes tracked with a bitset: includes all nodes in the graph but is updated locally
    bitset_t visited = bitset_create(n);

    // Local distance vector - stores shortest path distance for each vertex from source
    int *d = malloc(g.n_vertices * sizeof(int));
    for (int i = 0; i < g.n_vertices; i++)
        d[i] = INF;

    // Local parent array - stores parent of each vertex in BFS tree
    int *parent = malloc(g.n_vertices * sizeof(int));
    for (int i = 0; i < g.n_vertices; i++)
        parent[i] = -1;

    // Local frontier
    bitset_t frontier = bitset_create(n);
    node_t frontier_size = 0;


    start_time = MPI_Wtime();
    // Initialize: rank owning the source sets distance to 0 and adds to frontier
    if (source >= local_start && source < local_start + g.n_vertices) {
        node_t u = source - local_start;
        bitset_set(&visited, source);
        d[u] = 0;
        parent[u] = source;              // Source is its own parent
        bitset_set(&frontier, source);
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
        // Explore local frontier: for each frontier vertex, calculate send_counts for each rank
        for (node_t i = 0; i < local_start+g.n_vertices; i++) {
            if(bitset_test(&frontier, i)) {
                node_t u = i - local_start;                     // node in frontier - normalized
                for (node_t e = g.row_ptr[u] - local_edges; e < g.row_ptr[u + 1] - local_edges; e++) { //for all neighbors of u (normalized)
                    node_t v = g.col_idx[e];                    // neighbor of node in frontier (global)
                    if (bitset_test(&visited, v) == 1) {continue;}
                    int owner = v / block_size;                 // process owning the destination vertex of the edge
                    send_counts[owner]++;
                }
            }
        }
        
        // Allocate send_buf based on calculated send_counts -> avoiding over allocation of n_edges for all ranks...
        // still an over-estimation: we're counting more than one times the same vertex, but we're going to send only 
        // one edge for that outgoing vertex
        for (int i = 0; i < p; i++) {
            send_buf[i] = malloc(send_counts[i] * sizeof(Edge));
            send_counts[i] = 0; // Reset send_counts for reuse in the next loop 
        }

        // Populate send_buf with edges to be sent
        for (node_t i = local_start; i < local_start+g.n_vertices; i++) {
            if(bitset_test(&frontier, i)) {
                node_t u = i - local_start;                // node in frontier - normalized
                for (node_t e = g.row_ptr[u] - local_edges; e < g.row_ptr[u + 1] - local_edges; e++) {
                    node_t v = g.col_idx[e];            // neighbor of node in frontier (global)
                    if (bitset_test(&visited, v) == 1) {continue;}
                    bitset_set(&visited, v);

                    int owner = v / block_size;         // process owning the destination vertex of the edge
                    
                    send_buf[owner][send_counts[owner]].src = i;
                    send_buf[owner][send_counts[owner]++].dst = v;
                }
            }
        }

        /* ==============================================================
        Exchange of NUMBER of outgoing edges before Alltoallv -> all ranks know how much edges it needs to expect
        ============================================================== */
        int *recv_counts = malloc(p * sizeof(int));
        MPI_Alltoall(send_counts, 1, MPI_INT,
                     recv_counts, 1, MPI_INT,
                     comm);
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
        Exchange of ACTUAL outgoing edges before Alltoallv -> all ranks can now update their frontier
        ============================================================== */
        MPI_Alltoallv(send_flat, send_counts, sdispls, MPI_EDGE,
                      recv_buf, recv_counts, rdispls, MPI_EDGE,
                      comm);
    
        // Build next frontier: add unvisited vertices and set parent pointers
        frontier_size = 0;
        bitset_clear(&frontier);
        for (int i = 0; i < total_recv; i++) {
            Edge e = recv_buf[i];
            if (e.dst < local_start || e.dst > local_start + g.n_vertices) {
                printf("wtf: %d->%d arrived in %d\n", e.src, e.dst, rank);
                continue;
            }
            node_t local_v = e.dst - local_start;
            if (d[local_v] == INF) {
                d[local_v] = level;
                parent[local_v] = e.src;            // Set parent to the vertex that discovered this one
                bitset_set(&frontier, e.dst);    // Insert into next frontier
                frontier_size++;
                //printf("LEVEL %d - %d->%d arrived in %d, now having frontier size of %d\n", level, e.src, e.dst,rank, frontier_size);
            }
        }
        //printf("LEVEL %d - Rank %d with frontier of size %d\n", level, rank, frontier_size);

        /* ==============================================================
        Check global termination: stop if no new vertices were discovered
        ============================================================== */ 
        node_t global_frontier;
        MPI_Allreduce(&frontier_size, &global_frontier, 1,
                      MPI_UINT32_T, MPI_SUM, comm);

        if (global_frontier == 0) {
            end_time = MPI_Wtime();
            printf("rank%d spent %.3f seconds\n", rank, end_time-start_time);
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
            if(rank == 0) printf("BFS finished at level %d\n", level);            
            break;
        }

        level++;

    }


    /* 
    gathering of global parent array
    TODO: write this on a separate file with MPI-IO
    */
    int *recvcounts = NULL;
    int *displs = NULL;
    int *parent_global = NULL;

    if (rank == 0) {
        recvcounts = malloc(p * sizeof(int));
        displs = malloc(p * sizeof(int));

        // rank 0 can easily compute the offset of each parent array:
        for (int r = 0; r < p; r++) {
            int r_start = r * block_size;
            int r_end   = r_start + block_size;
            if (r_end > n) r_end = n;

            recvcounts[r] = (r_start < n) ? (r_end - r_start) : 0;
            displs[r] = r_start;

            parent_global = malloc(n * sizeof(int));
        }
    }

    MPI_Gatherv(parent, g.n_vertices, MPI_INT, parent_global, recvcounts, displs, MPI_INT, 0, comm);
    if (rank == 0) {
        printf("TOTAL PARENT ARRAY = ");
        int cap = 30 < n ? 30 : n;
        for (int i=0; i<cap; i++) {
            printf("%d ", parent_global[i]);
        } printf("\n");
    }

    free(send_buf);
    free(send_counts);
    bitset_free(&visited);
    bitset_free(&frontier);
    free(d);
    free(parent);

    return;
}

/* ============================================================================================================================
============================================================== MAIN ===========================================================
============================================================================================================================ */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    init_mpi_datatypes();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    CSR graph;
    graph.row_ptr = NULL;
    graph.col_idx = NULL;
    graph.n_vertices = 0;
    graph.n_edges = 0;

    // TODO: either randomize or selected as a program argument...
    node_t source = 0;
    edge_t start, end;
    node_t total_vertices;

    /* ============================================================== 
    STEP 1: each rank reads its portion of CSR through MPI IO
    ============================================================== */
    MPI_File handle;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &handle);

    // manually making the node_t and edge_t types match with MPI datatypes
    MPI_File_read_at(handle, 0, &total_vertices, 1, MPI_UINT32_T, NULL);

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

    start = graph.row_ptr[0];
    end = graph.row_ptr[graph.n_vertices];
    graph.n_edges = end-start;          //each process will own all the outgoing edges of the nodes they own
    graph.col_idx = malloc(graph.n_edges * sizeof(node_t));

    long long col_offset = sizeof(node_t) + sizeof(edge_t) + (total_vertices+1) * sizeof(edge_t) + (start) * sizeof(node_t);

    MPI_File_read_at(handle, col_offset, graph.col_idx,   
                graph.n_edges, MPI_UINT32_T, NULL);

    MPI_File_close(&handle);

    #ifdef DEBUG
    printf("%d with local start = %d\n", rank, rank*block_size);
    print_csr_stats(&graph);
    #endif

    printf("RANK %d - %d vertices and %ld edges, going from %ld to %ld\n", rank, graph.n_vertices, graph.n_edges, start, start + graph.n_edges);
    //print_csr_stats(&graph);


    /* ============================================================== 
    STEP 2: BFS algorithm
    ============================================================== */
    distributed_bfs(graph, rank*block_size, total_vertices, source, MPI_COMM_WORLD);


    //TODO: graph500 checks on parent array

    free_csr(&graph);

    MPI_Type_free(&MPI_EDGE);
    MPI_Finalize();
    return 0;
}