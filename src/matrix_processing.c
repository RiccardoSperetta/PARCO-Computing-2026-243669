#include <linux/limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

int ROWS, COLS, nnz;

// data will be stored directly in a triplet
typedef struct {
    int row;
    int col;
    double val;
} Triplet;

// -> easier ordering
// function to pass as an argument to qsort: how to order a matrix
int compare_by_row(const void *a, const void *b) {
    const Triplet *t1 = (const Triplet *)a;
    const Triplet *t2 = (const Triplet *)b;
    if (t1->row != t2->row)
        return (t1->row - t2->row);
    return (t1->col - t2->col); // secondary sort by column
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s matrix_name\n", argv[0]);
        return 1;
    }

    char* matrix_name = argv[1];

/*  ======= SETUP =============================================
    =========================================================== */
    //builds full file path:
    char input_path[256];
    snprintf(input_path, sizeof(input_path), "data/raw/%s.mtx", matrix_name);

    FILE *fp = fopen(input_path, "r");
    if (!fp) {
        perror("fopen");
        fprintf(stderr, "Failed to open: %s\n", input_path);
        return 1;
    }

    //arrays to be filled: full matrix and the 3 separate CSR arrays to compute
    Triplet* matrix;
    int* RowPtr;
    int* Acol;
    double* Aval;


/*  ======= EXTRACT MATRIX ====================================
    =========================================================== */
    //variables in order to read from file correctly
    char *line = NULL;
    size_t len = 0;
    int nread;
    bool matrix_data = true;
    bool matrix_header = true;
    int iteration = 0;
    bool symmetric = false;

    while ((nread = getline(&line, &len, fp)) != -1) {
        //parsing the matrix header:
        if (matrix_header && nread > 0) {
            //reads if matrix is real
            if (strcasestr(line, "real") == NULL) {
                fprintf(stderr, "Unsupported matrix type: only real matrices are supported.\n");
                free(line);
                fclose(fp);
                return 1;
            }
            //reads if matrix is symmetric
            symmetric = (strcasestr(line, "symmetric") != NULL);
            
            matrix_header = false;
            continue;
        }
        
        // skip comment lines starting with '%'
        if (nread > 0 && line[0] == '%') continue;
        
        // get the general data of the matrix
        if (matrix_data && nread > 0) {
            if (sscanf(line, "%d %d %d", &ROWS, &COLS, &nnz) == 3) {
                if (symmetric) {
                    nnz = nnz*2; //ideally will hold double of the entries mentioned in the first line
                }
                matrix = (Triplet*) malloc(nnz*sizeof(Triplet));
            } else {
                fprintf(stderr, "bad formatting of mtx file...\n");
                return 1;
            }

            matrix_data = false; //once found, look for the matrix contents:
        } 
        
        //fill the matrix: COO format
        else {
            Triplet t;
            sscanf(line, "%d %d %lf", &t.row, &t.col, &t.val);
            if (t.val == 0) {
                if (symmetric) { //there wont be no specular values stored in the matrix
                    nnz--;
                }
                nnz--; //some matrices also include "0" entries -> ignored 
            } else {
                //because indexes are 1 based:
                t.col--;
                //not doing t.row--: makes the CSR computation easier later
                matrix[iteration] = t;
                iteration++;
                if (symmetric) {
                    if (t.row == (t.col+1)) { //value is on the matrix diagonal
                        nnz--;                // = no specular will be stored          
                    } else {
                        Triplet specular = {t.col+1, t.row-1, t.val};
                        matrix[iteration] = specular;
                        iteration++;
                    }
                }
                
            }
        }
    }

    free(line);
    fclose(fp);


/*  ======= BUILDING CSR ======================================
    =========================================================== */
    //SORTING matrix by ROWS:
    qsort(matrix, nnz, sizeof(Triplet), compare_by_row);


    #ifdef DEBUG
    printf("COO matrix:\n");
    for (int i = 0; i<20; i++) {
        //even if stored correctly, if a number is too small, it wont be printed simply with %lf
        printf("%d %d %.15e\n", matrix[i].row - 1, matrix[i].col, matrix[i].val);
    }
    printf("...\n");    
    #endif

    RowPtr = (int*) calloc((ROWS+1), sizeof(int));
    Acol = (int*) malloc(nnz*sizeof(int));
    Aval = (double*) malloc(nnz*sizeof(double));

    // fill the array of pointers to the start of each row:
    int prevRow = 0;
    for(int i=0; i<nnz; i++) {
        int nextRow = matrix[i].row; //remember from above: here I have the row index + 1 (= the pointer for the start of the next row)
        if(nextRow>prevRow) {
            //if I'm jumping multiple rows I fill them with the previous value
            //= all those rows start and finish at the same index = contain 0 elements
            for(int j=prevRow+1; j<=nextRow; j++){
                RowPtr[j] = RowPtr[prevRow];
            }
        }
        RowPtr[nextRow]++;
        //fillign Acol and Aval is straight forward
        Acol[i] = matrix[i].col;
        Aval[i] = matrix[i].val;

        prevRow = nextRow;
    }

    free(matrix);


/*  ======= STORING CSR =======================================
    =========================================================== */
    //creating path for processed matrix data:
    char rowptr_path[PATH_MAX], col_path[PATH_MAX], val_path[PATH_MAX];
    snprintf(rowptr_path, sizeof(rowptr_path), "data/processed/%s/rowptr.bin", matrix_name);
    snprintf(col_path, sizeof(col_path), "data/processed/%s/col.bin", matrix_name);
    snprintf(val_path, sizeof(val_path), "data/processed/%s/val.bin", matrix_name);


    // opening in write mode the paths, 
    // will create the files if they don't exist already, and overwrite them otherwise:
    FILE *fp_rowptr = fopen(rowptr_path, "wb");
    FILE *fp_col = fopen(col_path, "wb");
    FILE *fp_val = fopen(val_path, "wb");

    //if something went wrong:
    if (!fp_rowptr || !fp_col || !fp_val) {
        perror("fopen output files");
        // Close any that succeeded
        if (fp_rowptr) fclose(fp_rowptr);
        if (fp_col) fclose(fp_col);
        if (fp_val) fclose(fp_val);
        return 1;
    }

    //writing down the data stored in arrays + matrix size values
    fwrite(&ROWS, sizeof(int), 1, fp_rowptr);
    fwrite(RowPtr, sizeof(int), ROWS + 1, fp_rowptr);

    fwrite(&COLS, sizeof(int), 1, fp_col);
    fwrite(Acol, sizeof(int), nnz, fp_col);

    fwrite(&nnz, sizeof(int), 1, fp_val);
    fwrite(Aval, sizeof(double), nnz, fp_val);

    fclose(fp_rowptr);
    fclose(fp_col);
    fclose(fp_val);


/*  ======= FINAL FREES =======================================
    =========================================================== */
    free(RowPtr);
    free(Acol);
    free(Aval);
}