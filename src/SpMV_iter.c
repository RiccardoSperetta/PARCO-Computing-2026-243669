#include <linux/limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "timer.h" //header with time measurement tools

#ifdef _OPENMP
#include <omp.h>
#endif

# ifdef DEBUG
void print_parallel_info(void); //debug print
# endif

// matrix sizes: rows, columns, number of non-zero elements
int ROWS, COLS, nnz;


int main(int argc, char* argv[]) {    
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <matrix_name>\n", argv[0]);
        return 1;
    }

    char *matrix_name = argv[1];

/*  ======= INITIALIZATION ====================================
    =========================================================== */
    // Build paths to the binary files, 
    // representing the 3 vectors that will be used in the multiplication based on CSR 
    char rowptr_path[PATH_MAX];
    char col_path[PATH_MAX];
    char val_path[PATH_MAX];

    snprintf(rowptr_path, sizeof(rowptr_path), "data/processed/%s/rowptr.bin", matrix_name);
    snprintf(col_path, sizeof(col_path), "data/processed/%s/col.bin", matrix_name);
    snprintf(val_path, sizeof(val_path), "data/processed/%s/val.bin", matrix_name);

    FILE *fp_rowptr = fopen(rowptr_path, "rb");
    FILE *fp_col = fopen(col_path, "rb");
    FILE *fp_val = fopen(val_path, "rb");

    if (!fp_rowptr || !fp_col || !fp_val) { //if anything goes wrong
        perror("fopen input files");
        return 1;
    }

    //initialization of sizes:
    //ASSUMPTION: no error handling since files are expected to be built with matrix_processing.c
    int trash;
    trash = fread(&ROWS, sizeof(int), 1, fp_rowptr);
    trash = fread(&COLS, sizeof(int), 1, fp_col);
    trash = fread(&nnz, sizeof(int), 1, fp_val);
    
    //RowPtr initialization:
    int *RowPtr = malloc((ROWS+1) * sizeof(int));
    if (!RowPtr) {
        fprintf(stderr, "malloc failed for RowPtr\n");
        fclose(fp_rowptr);
        return 1;
    }

    trash = fread(RowPtr, sizeof(int), ROWS+1, fp_rowptr);
    fclose(fp_rowptr);

    //Acol initialization:
    int* Acol = malloc(nnz * sizeof(int));
    if (!Acol) {
        fprintf(stderr, "malloc failed for Acol\n");
        fclose(fp_col);
        free(RowPtr);
        return 1;
    }

    trash = fread(Acol, sizeof(int), nnz, fp_col);
    fclose(fp_col);

    //Aval initialization:
    double* Aval = malloc(nnz * sizeof(double));
    if (!Aval) {
        fprintf(stderr, "malloc failed for Aval\n");
        fclose(fp_val);
        free(RowPtr);
        free(Acol);
        return 1;
    }

    trash = fread(Aval, sizeof(double), nnz, fp_val);
    fclose(fp_val);


    double now;
    GET_TIME(now);
    srand(now);     //simple randomization based on time

    //VECTOR to multiply the matrix with
    double* vector = malloc(COLS*sizeof(double));
    for (int i = 0; i<COLS; i++) {
        #ifdef DEBUG
            vector[i] = 1.; //more varied multiplications
        #else
            vector[i] = rand()/1000000.0 + rand()%100; //more varied multiplications
        #endif
    }

    //RESULT VECTOR of matrix x vector, all initialized with 0s
    double* result = calloc(COLS, sizeof(double));



/*  ======= MATRIX VECTOR MULTIPLICATION ======================
    =========================================================== */
    double start, finish;

    // multiple iterations in order to measure warm cache behavior
    for (int i = 0; i<11; i++) {
        //multiplication: either with sequential or parallel code:
        GET_TIME(start);
        #ifdef _OPENMP
#       pragma omp parallel for schedule(runtime) //schedule is determined at compile time
        #endif 
        for (int i=0; i<ROWS; i++) {
            double sum = 0.0;
            for(int j=RowPtr[i]; j<RowPtr[i+1]; j++){
                sum += Aval[j] * vector[Acol[j]];
            }
            result[i] = sum;
        }
        GET_TIME(finish);

        if(i>0) { // cold start is ignored
            double elapsed = finish-start;
            printf("%e\n", elapsed);
        }
    }

/*  ======= DEBUG PRINTS ======================================
    =========================================================== */
#ifdef DEBUG
    #ifdef _OPENMP
        printf("parallel code running...\n");
        print_parallel_info();
    #else
        printf("sequential code running...\n");
    #endif

    printf("rows= %d columns= %d nnz= %d\n", ROWS, COLS, nnz);
    //just for checking the first 10 elements
    int debug_print_size = 10;
    printf("RowPtr: ");
    for(int i=0; i<debug_print_size; i++) {
        printf("%d ", RowPtr[i]);
    }
    printf("\nAcol: ");
    for(int i=0; i<debug_print_size; i++) {
        printf("%d ", Acol[i]);
    }
    printf("\nAval: ");
    for(int i=0; i<debug_print_size; i++) {
        printf("%.2f ", Aval[i]);  //use %.15e for a more precise print
    }
    printf("\nVECTOR: ");
    for(int i=0; i<debug_print_size; i++) {
        printf("%lf ", vector[i]);
    }
    printf("\nRESULT: ");
    for(int i=0; i<debug_print_size; i++) {
        printf("%lf ", result[i]);
    }
    printf("\n");
#endif


/*  ======= FINAL FREES =======================================
    =========================================================== */
    free(RowPtr);
    free(Acol);
    free(Aval);

    free(vector);
    free(result);

    return 0;
}


/*  ======= end of main =======================================
    =========================================================== */


// debugging function: prints the schedule and chunksize being utilized:
#ifdef _OPENMP
void print_parallel_info(void) {
    omp_sched_t kind;
    int chunk_size;
    omp_get_schedule(&kind, &chunk_size);
    
    //printf("%d\n", kind);
    char* modifier;
    if((kind & 0x80000000) != 0) {
        modifier = "monotonic";
        kind = kind & 0x7FFFFFFF;  // Mask off the monotonic bit
    } else {
        modifier = "nonmonotonic";
    }

    // Match it
    const char *kind_str;
    if (kind == omp_sched_static) {
        kind_str = "static";
    } else if (kind == omp_sched_dynamic) {
        kind_str = "dynamic";
    } else if (kind == omp_sched_guided) {
        kind_str = "guided";
    } else if (kind == omp_sched_auto) {
        kind_str = "auto";
    } else {
        kind_str = "UNKNOWN";
    }
    
    printf("CORES: %d; THREADS: %d; SCHEDULE: %s:%s; ", omp_get_num_procs(), omp_get_max_threads(), modifier, kind_str);
    if (chunk_size > 0)
        printf("CHUNKSIZE: %d\n", chunk_size);
    else
        printf("default chunk size (implementation-defined)\n");
}
#endif