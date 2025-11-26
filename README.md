## Table of Contents

1. [Project Overview](#1-project-overview)

2. [Repo layout](#2-repo-layout)

3. [Compilation](#3-compilation)

4. [Running the Code](#4-running-the-code)

5. [Input Data](#5-input-data)

6. [Modifying Parameters](#6-modifying-parameters)

7. [Cluster Notes](#7-cluster-notes)

8. [Reproducibility and Data](#8-reproducibility-and-data)

  
## 1. Project Overview

Author: Riccardo, Speretta 243669 riccardo.speretta@studenti.unitn.it
This project constitutes all the necessary tools to reproduce the data used in the Deliverable 1 for the Parallel Computing 2025/2026 course at unitn.

---
## 2. Repo layout

```bash

├── README.md

├── src/        # C/C++ source code

├── scripts/    # run/plot scripts

├── results/    # output logs/timing results

├── plots/      # figures for the report

├── paper_data/ # all utilized data produced by the cluster + relative graphs

└── data/       # input matrices

```

#### src
- `matrix_processing.c`
  = taken the name of a matrix in input this program will look for it in the `data/raw` folder(as `matrix_name.mtx`, and if it's present it will translate it into the CSR format, specifically it will create 3 binary files inside of `data/processed/matrix_name`:
    1) `rowptr.bin` = (1 long) = number of rows (ROWS), followed by the list of pointers to the start of each row in the columns array (ROWS+1 longs)

    2) `col.bin` = (1 long) = number of columns(COLS), followed by the list of column indexes (nnz longs)

    3) `val.bin` = (1 long) = nnz(non zero values), followed by the list of values (nnz doubles) contained in the matrix

- `SpMV_iter.c`
  = taken the name of a matrix in input this program will look for it in the `data/processed` folder, expecting the previous program to be run already.
  It will then perform the standard matrix/vector multiplication with CSR format, including the parallel openMP `#pragma` to allow for parallelization when compiled with `-fopenmp`
  The multiplication is performed 11 times in a row: the first run is to warm up the cache, the following are actually registered as values and printed

(*both of the previous programs include debugging sections, where extra prints are provided when compiled with the `-DDEBUG` flag*)

- `timer.h`
  = defines a simple function for time measurements
  based on the POSIX `gettimeofday()` function, which provides wall-clock time with microsecond resolution. The timer records the elapsed real time (including OS overhead and thread scheduling effects). Time values are stored as doubles in milliseconds, computed as `t.tv_sec + t.tv_usec / 1e6`. This ensures consistent and comparable timing across all experiments

#### scripts

- `matrix_list.sh`
    = the list of all the matrices that we intend using
    it's shared among the other scripts, so during a job submission it's best not to modify this file at all

- `download_matrix.sh`
    = given the name and group of the matrix, this script will download and unzip the corresponding matrix from https://sparse.tamu.edu/, storing it in `data/raw`

- `download_list.sh`
    = calls the previous script for all matrices in the list

- `pre_process_matrix.sh`
    = given the name of the matrix, calls `matrix_processing.c` on the downloaded `.mtx` matrix file, fails if file is missing.

- `pre_process_list.sh`
    = calls the previous script for all matrices in the list

- `measure_matrix.sh`
    = given the name of a matrix, performs the actual benchmarking with `SpMV_iter.c`, changing all environment variables iteratively for all possible configurations (of optimization levels for sequential code and thread count and omp schedules), performing 10 execution for time measurements, and 10 other separate measurements for cache performacne. 

- `measure_list.sh`
    = calls the previous script for all matrices in the list

- `batch.pbs`
    = the pbs script intended for a cluster job submission
    -> includes the execution of all the previous `*_list.sh` scripts in order to build an easy to understand and modular pipeline for benchmarking
    (if matrices are already present or processed, they wont be downloaded/processed again)

each script will provide helpful prints to allow the user understand the phase reached by the scripts in their execution.
  
#### scripts/plotting
nested folder dedicated to python scripts responsible for plotting:
- `loaders.py` 
  = reads data from the `results/` folder and produces data frames that python can work upon
- `summaries.py`
  = creates a summary with mean + 90th percentile for the given value columns
- `plot_results.py`
  = actually calls different functions to create the desired graphs and stores them in `plots/`

*note*: `requirements_py.txt`(in root folder) includes all specific dependencies used in the python environment, can be used as reference for which libraries are required to obtain identical graphs.
  
#### results
(its contents, given their sizes, are not included by default in the repo, but will be populated by the previous scripts)

for each matrix, like for *data/processed*, there will be one directory, containing:
- **sequential** = all optimization levels corresponding folders 

- **thx** = where x is a power of 2, indicating the number of threads utilized, 
  inside of each there will be all the `schedule_chunk-size` combinations between 
  (static, dynamic, guided) x (10, 100, 1000)

each "end point folder" contains `perf.txt` and `time.txt`, respectively for perf cache misses percentages (for LL1 and LLC) and time measurements.

#### plots
(meant to be filled locally with the Python scripts, not part of the cluster execution of the batch job)

#### paper_data
- `results.zip` include the original `results/` folder (zipped) that I've used for graph creation
- `plots/` is exactly the copy with all the pictures of all graphs (including un-utilized ones) of the paper
- `job_outputs` contains the `txt` files produced by the job submissions I've made to the cluster

#### data
(its contents, given their sizes, are not included by default in the repo, but will be populated by the previous scripts)
- **raw**
    = all the `matrix_name.mtx` files downloaded

- **processed**
    = contains the outputs of `matrix_processing.c`:
    the 3 binary files are stored in one directory per matrix( = matrix name without the .mtx),

  ---
## 3. Compilation
The compiler version used on the cluster for all compilations of C code is `gcc9.1.0`

As mentioned above all of the compilation is done automatically by the scripts, so the user has very little room for mistakes:
The `batch.pbs` loads the `gcc91` module, so that each script can override the standard `gcc` command with `gcc() { gcc-9.1.0 "$@"; }` and use the correct compiler version.

- `matrix_processing.c` is compiled with `-O3`
- `SpMV_iter.c` varies in its compilation accordingly to the iterative configuration built by the `measure_matrix.sh` script (a set of basic flags is shared between sequential and parallel execution, parallel code includes always `-O3` and obviously `-fopenmp`)

*plus*: every compilation requires `--std=c11` to work properly on the cluster environment 

---
## 4. Running the code
ALL scripts and executables are expected to be run inside the root folder of the project. 

The pipeline execution is reported in the `batch.pbs` script of the project:
- `download_list.sh` assures all matrices in the `matrix_list.sh` are downloaded correctly
  if no errors arise ->
- `pre_process_list.sh` will process those matrices
  if no errors arise ->
- `measure_list.sh` will actually perform all of the benchmarking for all matrices

The job can simply be submitted via `qsub batch.pbs`
At the same time, single `*_matrix.sh` scripts can be run to perform individual steps separately.

---
## 5. Input Data
Input data, as described above, can be easily modified:
- in terms of matrices those can be selected by reporting their name and group (must be valid and present in the [SuiteSparse matrix collection](https://sparse.tamu.edu/)) inside of `matrix_list.sh`
- clearly, those matrices can also be download manually, if correctly named and positioned inside of `data/raw`


---
## 6. Modifying Parameters
- number of CPUs and cluster related parameters can be changed from `batch.pbs` script, overriding the already existing ones, but though the results gathering should execute correctly the plotting of those data is not guaranteed to be working.
- number of threads is modified iteratively by `measure_matrix.sh` via
  `export OMP_NUM_THREADS=$threads`
- schedule and chunk-size is modified iteratively by `measure_matrix.sh` via
  `export OMP_SCHEDULE="${sched},${chunk}"`

---
## 7. Cluster Notes
Cluster-specific notes (modules, queues)
- be sure to be running the job on the `hpc-c11` nodes of the unitn cluster in order to obtain comparable performances 
- the utilized queue is `short_cpuQ`, and further specifications are included inside the script itself
  
---
## 8. Reproducibility and Data
- For full transparency the exact results, plots, and job outputs are also reported in the `paper_data` folder, this indeed includes all the data that was used for the data shown inside the paper + graphs that are not included in it for spatial reasons.
  This way reproducibility of the results obtained can be compared to the one I've used.

