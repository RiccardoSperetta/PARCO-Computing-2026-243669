#include <linux/limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "timer.h"

long ROWS, COLS, nnz;


int main(int argc, char* argv[]) {    
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <matrix_name>\n", argv[0]);
        return 1;
    }

    char *matrix_name = argv[1];

// INIT - TODO: move into init.h
    // Build paths to the binary files
    char rowptr_path[PATH_MAX];
    char col_path[PATH_MAX];
    char val_path[PATH_MAX];

    snprintf(rowptr_path, sizeof(rowptr_path), "data/processed/%s/rowptr.bin", matrix_name);
    snprintf(col_path, sizeof(col_path), "data/processed/%s/col.bin", matrix_name);
    snprintf(val_path, sizeof(val_path), "data/processed/%s/val.bin", matrix_name);

    FILE *fp_rowptr = fopen(rowptr_path, "rb");
    FILE *fp_col = fopen(col_path, "rb");
    FILE *fp_val = fopen(val_path, "rb");

    if (!fp_rowptr || !fp_col || !fp_val) {
        perror("fopen input files");
        return 1;
    }

    //initialization of sizes:
    fread(&ROWS, sizeof(long), 1, fp_rowptr);
    fread(&COLS, sizeof(long), 1, fp_col);
    fread(&nnz, sizeof(long), 1, fp_val);
#ifdef DEBUG
    printf("%ld %ld %ld\n", ROWS, COLS, nnz);
#endif

    //RowPtr initialization:
    long *RowPtr = malloc((ROWS+1) * sizeof(long));
    if (!RowPtr) {
        fprintf(stderr, "malloc failed for RowPtr\n");
        fclose(fp_rowptr);
        return 1;
    }

    fread(RowPtr, sizeof(long), ROWS+1, fp_rowptr);
    fclose(fp_rowptr);

    //Acol initialization:
    long* Acol = malloc(nnz * sizeof(long));
    if (!Acol) {
        fprintf(stderr, "malloc failed for Acol\n");
        fclose(fp_col);
        free(RowPtr);
        return 1;
    }

    fread(Acol, sizeof(long), nnz, fp_col);
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

    fread(Aval, sizeof(double), nnz, fp_val);
    fclose(fp_val);


    //vector(to multiply the matrix with) and result vector initialization:
    double* vector = malloc(COLS*sizeof(double));
    for (long i = 0; i<COLS; i++) {
        vector[i] = 1.; //TODO: randomize!
    }
    double* result = calloc(COLS, sizeof(double));


//SpMV WITH CSR ================================
    double start, finish;

    //multiplication: either with sequential or parallel code:
#ifdef _OPENMP
    GET_TIME(start);
#   pragma omp parallel for
    for (long i=0; i<ROWS; i++) {
        for(int j=RowPtr[i]; j<RowPtr[i+1]; j++){
            result[i] += Aval[j] * vector[Acol[j]];
        }
    }
    GET_TIME(finish);
#else 
    GET_TIME(start);
    for (long i=0; i<ROWS; i++) {
        for(int j=RowPtr[i]; j<RowPtr[i+1]; j++){
            result[i] += Aval[j] * vector[Acol[j]];
        }
    }
    GET_TIME(finish);
#endif
    double elapsed = finish-start;
    printf("%e\n", elapsed);

// DEBUG PRINT ==========================
#ifdef DEBUG
    #ifdef _OPENMP
        printf("parallel code running...\n");
    #else
        printf("sequential code running...\n");
    #endif

    printf("Elapsed time = %e seconds\n", elapsed);

    //just for checking the first 10 elements
    for(long i=0; i<10; i++) {
        printf("%lf ", result[i]);
    }
    printf("\n");
#endif

// FINAL FREES ==========================
    free(RowPtr);
    free(Acol);
    free(Aval);

    free(vector);
    free(result);

    return 0;
}