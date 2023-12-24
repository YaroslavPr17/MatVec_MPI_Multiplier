# Laboratory work #2. Linear algebra
**Implementation of the matrix-vector multiplication algorithm using the MPI library (Message Passing Interface)**

## Setting the task
Task: _implement 3 algorithms for multiplying a matrix by a vector_:
1. When splitting by lines
<p align="center">
<img src="https://github.com/YaroslavPr17/MatVec_MPI_Multiplier/assets/77925460/0b3ac7fd-de71-4d5e-b023-5f14b12721df" alt="image" width="50%" height="auto">
</p>

2. When divided into columns
<p align="center">
<img src="https://github.com/YaroslavPr17/MatVec_MPI_Multiplier/assets/77925460/f35792fb-4dff-4c9b-9819-a40afab27364" alt="image" width="70%" height="auto">
</p>

3. With block partitioning
<p align="center">
<img src="https://github.com/YaroslavPr17/MatVec_MPI_Multiplier/assets/77925460/cc371ccf-88b3-4ea3-948f-aa2f1fdc1bed" alt="image" width="40%" height="auto">
</p>


## Description of the source data
### Number of executing processes
* **The number of executing processes** is less than or a multiple of the number of physical CPU cores.
### Data type
* `double`
### Matrix
* **The size of the matrix** for consistency of experiments in each dimension is a multiple of each value of the number of executing processes
### Vector
* **The size of the vector** is equal to the second dimension of the matrix (the number of columns of the matrix). That is, if the size of the matrix is $(24*48)$, then the length of the vector will be $48$.

> Data generation was performed using the Python library `numpy`, followed by saving the values of `double` in the format `%.4f` in terms of C language specifiers to a text file.

> The matrices for performing calculations are square in order to preserve the total amount of data for each Process of a distributed system.

## Description of the output data
As a result of the algorithms, a vector of type `double` of the size of the first dimension of the matrix (the number of its rows) is expected, that is, $24$ under the condition above.

## Description of metrics
The following metrics were used to evaluate the performance of the algorithms
### Online
* **Execution time**. Conditions:
  * At the beginning of the operation time measurement, the algorithm has preloaded Matrix and Vector data on the main process.
  * Completion of the work time measurement occurs when the main process has received the final vector $\textemdash$ The result of multiplying the matrix by the vector.
  * The running time was measured independently on each process and was subsequently aggregated using the `max` function.
### Offline
* **Speed Up** (Acceleration)
  * It is calculated as a quotient of the execution time on one process and on `n` processes. $$S = \frac{T_{serial}}{T_{parallel}}$$
* **Efficiency** (Efficiency)
  * Let `S`$\textemdash$ be the acceleration, and `p` $\textemdash$ be the number of processes. Then the efficiency is calculated as a quotient of the acceleration and the number of processes. $$E = \frac{S}{p} = \frac{T_{serial}}{p * T_{parallel}}$$

> Averaging of 100 experiments was used for measurement accuracy in each experiment

## Computer configuration
* **CPU**: Intel(R) Core(TM) i5-10400F CPU @ 2.90GHz 2.90 GHz. | Cores: 6 | Hyper Threading
* **RAM**: 32.0 GB
* **OS**: WSL Ubuntu 22.04.3 LTS

## Visualization of the received metrics
### For each Algorithm (Line by line)
<p align="center">
<img src="https://github.com/YaroslavPr17/MatVec_MPI_Multiplier/assets/77925460/a186cabb-9c47-4d81-bc41-64148c8824db" alt="image" width="100%" height="auto">
</p>

### Algorithm comparison
<p align="center">
<img src="https://github.com/YaroslavPr17/MatVec_MPI_Multiplier/assets/77925460/1590242f-0a23-4e45-a009-4402630b2dcd" alt="image" width="100%" height="auto">
</p>


## Conclusions
### Split by rows
* The more data there is, the more noticeable the decrease in **work time**
* Processes up to the 6th (for 6 physical processor cores) are "useful". With further increase (up to 12) ** execution time increases moderately (for small amounts of data) and increases dramatically for any amount of data when increasing processes to 24. The reason for $\textemdash$ is the overhead of interprocess communication.
* **Acceleration** maximum for 6 processes (for 6 physical cores).
* The worst **acceleration** is when there is extremely little data for each process.
* **Acceleration** behaves the same for a larger proportion of datasets.
* **Efficiency** decreases exponentially depending on the number of processes.
* The least **efficiency** $\textemdash$ when dividing a small task into many Workers.

### Split by column (Differences from the row-split algorithm)
* The minimum **execution time** corresponds to 6 processes (for 6 physical CPU cores).
* The minimum data set has one of the best **speedups** for all processes.

### Split by blocks (Differences from the Algorithm with row-split)
* The best result **of time** on 12 processes (for 6 physical cores) with a large amount of data. The multithreading of the processor turns out to be useful.

### Algorithm comparison
* The value of **efficiency** and **acceleration** is better when divided into Columns under all conditions, however **execution time** is much worse than that of other algorithms.
* The shortest **time** is achieved using a Block Division algorithm.
* When divided into Columns depending on the amount of data **acceleration** and **efficiency** is 2-4 times better with 24 processes and 2 times better with 6 active processes (by 6 physical processor cores).
