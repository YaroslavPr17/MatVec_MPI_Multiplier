#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "matr_utils.h"
#include "utils.h"
#include "constants.h"



void distribute_data (double* matrix, double* vector, long int n_rows, long int n_cols, long int local_n, int my_rank, int comm_sz, MPI_Comm comm, double* local_matr){
    int error;

    if (my_rank == MAIN_PROCESS){
        error = MPI_Scatter(
            &matrix[0],
            local_n * n_cols,
            MPI_DOUBLE,
            local_matr,
            local_n * n_cols,
            MPI_DOUBLE,
            MAIN_PROCESS,
            comm
        );
    }
    else {
        error = MPI_Scatter(
            NULL,
            -1,
            MPI_DOUBLE,
            local_matr,
            local_n * n_cols,
            MPI_DOUBLE,
            MAIN_PROCESS,
            comm
        );
    }
    process_error(error);

    error = MPI_Bcast(
        vector,
        n_cols,
        MPI_DOUBLE,
        MAIN_PROCESS,
        comm
    );
    process_error(error);

    return;
}


int main(int argc, char** argv){
    int error;

    long int n_rows = strtol(argv[1], NULL, 10);
    long int n_cols = strtol(argv[2], NULL, 10);

    int comm_sz;
    int my_rank;

    double start, finish, elapsed = 0;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == MAIN_PROCESS){
        if (n_rows % comm_sz != 0){
            printf("\nERROR!!!\n%ld mod %d = %ld. Unable to parallellize task.\n", n_rows, comm_sz, n_rows % comm_sz);
            return 0;
        }

        char* new_file_name = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
        sprintf(new_file_name, "rowwise_%ld_%ld.csv", n_rows, n_cols);

        if (fopen(new_file_name, "r") == NULL){
            FILE* fp = fopen(new_file_name, "w");
            if (fp == NULL){
                printf("Unable to create output file.\n");
                return 0;
            }
            fprint(fp, "n_rows, n_cols, n_processes, time\n");
        }
    }

    long int local_n = n_rows / comm_sz;

    double* matrix;
    double* result;
    double* vector = (double*) malloc(n_cols * sizeof(double));


    if (my_rank == MAIN_PROCESS){
        printf("n_rows = %ld\n", n_rows);
        printf("n_cols = %ld\n", n_cols);
        printf("local_n = %ld\n", local_n);
        printf("comm_sz = %d\n", comm_sz);
        printf("my_rank = %d\n", my_rank);

        matrix = (double*) malloc(n_rows * n_cols * sizeof(double));
        result = (double*) malloc(n_rows * sizeof(double));


        if (my_rank == MAIN_PROCESS){
            error = load_matr(n_rows, n_cols, matrix);
            if (error == -1){
                char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
                build_matrix_filename(n_rows, n_cols, filename);
                printf("Unable to locate matrix file '%s'\n", filename);
                free(filename);
                return 0;
            }
            
            // print_matr(matrix, n_rows, n_cols, -1);
        }

        if (my_rank == MAIN_PROCESS){
            error = load_vec(n_cols, vector);
            if (error == -1){
                char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
                build_vector_filename(n_cols, filename);
                printf("Unable to locate vector file '%s'\n", filename);
                free(filename);
                return 0;
            }
            
            // print_vec(vector, n_cols, -1);
        }
    }
    // if (my_rank == MAIN_PROCESS){
    //     print_matr(matrix, n_rows, n_cols, my_rank);
    //     print_vec(vector, n_cols, my_rank);
    // }

    double* local_matr = (double*) malloc(local_n * n_cols * sizeof(double));
    double* local_result = (double*) malloc(local_n * sizeof(double));

    distribute_data(matrix, vector, n_rows, n_cols, local_n, my_rank, comm_sz, MPI_COMM_WORLD, local_matr);
    // print_matr(local_matr, local_n, n_cols, my_rank);
    // print_vec(vector, n_cols, my_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    multiply_std_rowwise(local_matr, vector, local_n, n_cols, local_result);
    // print_matr(local_matr, local_n, n_cols, my_rank);
    // print_vec(local_result, local_n, my_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    double local_elapsed = finish - start;
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, MAIN_PROCESS, MPI_COMM_WORLD);

    MPI_Gather(local_result, local_n, MPI_DOUBLE, result, local_n, MPI_DOUBLE, MAIN_PROCESS, MPI_COMM_WORLD);

    free(vector);
    free(local_matr);
    free(local_result);

    MPI_Finalize();

    if (my_rank == MAIN_PROCESS){
        print_vec(result, n_rows, -1);
        // printf("Elapsed time: %f s\n", elapsed);
        char* new_filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
        sprintf(new_filename, "rowwise_%ld_%ld.csv", n_rows, n_cols);

        FILE* fp = fopen(new_filename, "a+");
        if (fp == NULL){
            printf("Unable to open output file.\n");
            return 0;
        }
        fprint(fp, "n_rows, n_cols, n_processes, time\n");
        fprintf(fp, "%ld, %ld, %d, %lf\n");

        free(matrix);
        free(result);
    }

    return 0;
}

