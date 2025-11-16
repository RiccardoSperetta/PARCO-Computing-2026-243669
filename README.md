## Table of Contents
1. [Project Overview](#1-project-overview)
2. [Compilation](#2-compilation)
3. [Running the Code](#3-running-the-code)
4. [Input and Output](#4-input-data)
5. [Modifying Parameters](#5-modifying-parameters)
6. [Cluster Notes](#6-cluster-notes)
7. [Reproducibility and Data](#7-reproducibility-and-data)
8. [Results and Plots](#8-results-and-plots)

## 1. Project Overview
Author: Riccardo, Speretta 243669 riccardo.speretta@studenti.unitn.it




## 2. Repo layout
```bash
├── README.md
├── src/        # C/C++ source code
├── scripts/    # run/plot scripts
├── results/    # output logs/timing results
├── plots/      # figures for the report
└── data/       # input matrices 
```

#### src
- `matrix_processing.c`
= taken the name of a matrix in input will look for it in the data/raw folder, and if it's present it will translate it into the CSR format: ...
- `SpMV.c`
= taken the name of a matrix in input will look for it in the data/processed folder, 
and if it's present in the correct format it will load the 3 separate arrays representing the CSR format. It will then proceed to perform the matrix vector multiplication
(contains both the sequential and parallel code, since they can be differentiated by the -fopenmp flag at compile time)
- `timer.h`
= contains the basic method to measure time used inside of SpMV.c

#### scripts
- `matrix_list.sh`
    = the list of all the matrices that we intend using
    it's shared among the other scripts
- `download_matrix.sh`
    = given the name and group of the matrix, this script will download the corresponding matrix from https://sparse.tamu.edu/
- `download_list.sh`
    = calls the previous script for all matrices in the list
- `pre_process_matrix.sh`
    = s
- `pre_process_list.sh`
    = calls the previous script for all matrices in the list
- `measure_matrix.sh`
    = s
- `measure_list.sh`
    = calls the previous script for all matrices in the list

.
- `batch.pbs`
    = the pbs script intended for a cluster batch submission
    -> includes the execution of all the previous *_list.sh scripts in order to build an easy to understand and modular pipeline


#### data
(its contents, given their sizes, are not included by default in the repo, but will be populated by the previous scripts)
- **raw**
    = all the .mtx files downloaded
- **processed**
    = contains the outputs of matrix_processing.c:
    one directory per matrix( = matrix name without the .mtx), 
    here 3 binary files will be stored representing
    1) rowptr.bin = (1 long) = number of rows (ROWS), followed by the list of pointers to the start of each row in the columns array (ROWS+1 longs)
    2) col.bin = (1 long) = number of columns(COLS), followed by the list of column indexes (nnz longs)
    3) val.bin = (1 long) = nnz(non zero values), followed by the list of values (nnz doubles) contained in the matrix 



#### results
(its contents, given their sizes, are not included by default in the repo, but will be populated by the previous scripts)

for each matrix, like for *data/processed*, there will be one directory, containing:
- **sequential** = 
- **thx** = 

inside of each of those final directories there will be 2 files: 
1) time.txt = 
2) perf.txt = 


#### plots
(meant to be filled locally by the Python script, not part of the cluster execution of the batch job)

reading from the previous results, the Python script will compute different plots:




## 2. Compilation
TODO


## 3. Running the code
TODO



## 4. Input Data
This project uses the following matrices from the SuiteSparse Matrix Collection:

- still to decide definitely

They are downloaded from: https://sparse.tamu.edu/
(either manually or directly with download_matrix.sh script)




## 5. Modifying Parameters
what parameters actually?


## 6. Cluster Notes
TODO


## 7. Reproducibility and Data
The modularity of the repo should allow for an easy automation in all the steps of the pipeline,

for further details: it's possible to add the flag `-DDEBUG` in the compilation of SpMV.c in order to have a better understanding of what is being executed, 
also useful to check for the first 10 ... 


## 8. Results and Plots
TODO