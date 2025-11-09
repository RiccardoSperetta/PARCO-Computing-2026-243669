#include <linux/limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

long ROWS, COLS, nnz;

typedef struct {
    long row;
    long col;
    double val;
} Triplet;

int compare_by_row(const void *a, const void *b) {
    const Triplet *t1 = (const Triplet *)a;
    const Triplet *t2 = (const Triplet *)b;
    if (t1->row != t2->row)
        return (t1->row - t2->row);
    return (t1->col - t2->col); // secondary sort by column
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s file.mtx\n", argv[0]);
        return 1;
    }

    //builds full file path:
    char input_path[256];
    snprintf(input_path, sizeof(input_path), "data/raw/%s", argv[1]);

    FILE *fp = fopen(input_path, "r");
    if (!fp) {
        perror("fopen");
        fprintf(stderr, "Failed to open: %s\n", input_path);
        return 1;
    }

    char *line = NULL;
    size_t len = 0;
    long nread;

    Triplet* matrix;
    long* RowPtr;
    long* Acol;
    double* Aval;

    bool matrix_data = true;
    long iteration = 0;

//EXTRACT MATRIX ================================
    while ((nread = getline(&line, &len, fp)) != -1) {
        /* skip comment lines starting with '%' */
        if (nread > 0 && line[0] == '%') continue;
        /* get the general data of the matrix */
        if (matrix_data) {
            if (sscanf(line, "%ld %ld %ld", &ROWS, &COLS, &nnz) == 3) {
                matrix = (Triplet*) malloc(nnz*sizeof(Triplet));
                RowPtr = (long*) malloc((ROWS+1)*sizeof(long));
                Acol = (long*) malloc(nnz*sizeof(long));
                Aval = (double*) malloc(nnz*sizeof(double));
            } else {
                fprintf(stderr, "bad formatting of mtx file...\n");
                return 1;
            }

            matrix_data = false;
        } 
        //fill the matrix
        else {
            Triplet t;
            sscanf(line, "%ld %ld %lf", &t.row, &t.col, &t.val);
            //because indexes are 1 based:
            t.col--;
            //not doing t.row-- for simplicity in CSR computation
            matrix[iteration] = t;
            iteration++;
        }
    }

    free(line);
    fclose(fp);

//BUILDING CSR ================================
    //SORTING matrix by ROWS:
    qsort(matrix, nnz, sizeof(Triplet), compare_by_row);

    // fill the array of pointers to the start of each row:
    long prevRow = 0;
    for(long i=0; i<nnz; i++) {
        long nextRow = matrix[i].row;
        if(nextRow>prevRow) {
            for(long j=prevRow+1; j<=nextRow; j++){
                RowPtr[j] = RowPtr[prevRow];
            }
        }
        RowPtr[nextRow]++;
        Acol[i] = matrix[i].col;
        Aval[i] = matrix[i].val;

        prevRow = nextRow;
    }

    free(matrix);

//STORING PROCESSED MATRIX ================================
    //removing the .mtx:
    char basename[256];
    snprintf(basename, sizeof(basename), "%s", argv[1]);
    char *dot = strrchr(basename, '.');
    if (dot) *dot = '\0'; 

    //creating path for processed matrix data:
    char rowptr_path[PATH_MAX], col_path[PATH_MAX], val_path[PATH_MAX];
    snprintf(rowptr_path, sizeof(rowptr_path), "data/processed/%s/rowptr.bin", basename);
    snprintf(col_path, sizeof(col_path), "data/processed/%s/col.bin", basename);
    snprintf(val_path, sizeof(val_path), "data/processed/%s/val.bin", basename);

    //opening in write mode the paths:
    FILE *fp_rowptr = fopen(rowptr_path, "wb");
    FILE *fp_col = fopen(col_path, "wb");
    FILE *fp_val = fopen(val_path, "wb");

    if (!fp_rowptr || !fp_col || !fp_val) {
        perror("fopen output files");
        return 1;
    }

    //writing down the data stored in arrays + matrix size values
    fwrite(&ROWS, sizeof(long), 1, fp_rowptr);
    fwrite(RowPtr, sizeof(long), ROWS + 1, fp_rowptr);

    fwrite(&COLS, sizeof(long), 1, fp_col);
    fwrite(Acol, sizeof(long), nnz, fp_col);

    fwrite(&nnz, sizeof(long), 1, fp_val);
    fwrite(Aval, sizeof(double), nnz, fp_val);

    fclose(fp_rowptr);
    fclose(fp_col);
    fclose(fp_val);

//FINAL FREES ================================
    free(RowPtr);
    free(Acol);
    free(Aval);
}