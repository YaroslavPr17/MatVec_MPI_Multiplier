#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LENGTH 128
#define MAIN_PROCESS 0


char* build_matrix_filename(long int n_rows, long int n_cols){
    char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    sprintf(filename, "matrix_%ld_%ld.txt", n_rows, n_cols);

    return filename;
}

char* build_vector_filename(long int n_elems){
    char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    sprintf(filename, "vector_%ld.txt", n_elems);

    return filename;
}


void print_vec(double* vec, long int n_elems, int proc_num){
    printf("%d\tVECTOR:\n", proc_num);
    for (long int i = 0; i < n_elems; ++i){
        printf("%lf\n", vec[i]);
    }
    return;
}


int read_vector_from_file(long int n_elems, int my_rank, int comm_sz, MPI_Comm comm, double* vector){

    if (my_rank == MAIN_PROCESS){
        char* body = build_vector_filename(n_elems);
        char filename[MAX_FILENAME_LENGTH] = "./data/";
        strcat(filename, body);
        free(body);
        printf("Reading vector from file '%s'...\n", filename);

        FILE *fp = fopen(filename, "r");
        if (fp == NULL){ 
            return -1;
        }

        for (long int i = 0; i < n_elems; ++i){
            fscanf(fp, "%lf", &vector[i]); 
        }
    }

    int error;

    // printf("%d, %d, %d, %d, %d, %d\n", 
    //     MPI_SUCCESS, MPI_ERR_BUFFER, MPI_ERR_COMM,
    //     MPI_ERR_COUNT, MPI_ERR_TYPE, MPI_ERR_OTHER);

    MPI_Bcast(
        vector,
        n_elems,
        MPI_DOUBLE,
        MAIN_PROCESS,
        comm
    );

    // if (error != MPI_SUCCESS){
    //     printf("Error %d\n", error);

    //     char* str = malloc(128 * sizeof(char));
    //     int str_len; 

    //     MPI_Error_string(error, str, &str_len);

    //     printf("Error: %s, length: %d\n", str, str_len);

    // }

    return 0;
}


void print_matr(double* matr, long int n_rows, long int n_cols, int proc_num){
    printf("%d\tMATRIX:\n", proc_num);
    for (long int i = 0; i < n_rows; ++i){
        for (long int j = 0; j < n_cols; ++j){
            printf("%lf ", matr[i * n_cols + j]);
        }   
        printf("\n");
    }
    return;
}


int read_matrix_from_file(double* local_matr, long int local_n, long int n_rows, long int n_cols, int my_rank, int comm_sz, MPI_Comm comm){
    double* matrix;

    if (my_rank == MAIN_PROCESS){
        char* body = build_matrix_filename(n_rows, n_cols);
        char filename[MAX_FILENAME_LENGTH] = "./data/";
        strcat(filename, body);
        free(body);
        printf("Reading matrix from file '%s'...\n", filename);

        matrix = (double*) malloc(n_rows * n_cols * sizeof(double));

        FILE *fp = fopen(filename, "r");
        if (fp == NULL){
            return -1;
        }

        for (long int i = 0; i < n_rows; ++i){
            for (long int j = 0; j < n_cols; ++j){
                fscanf(fp, "%lf", &matrix[i * n_cols + j]);
            }   
        }

        // print_matr(matrix, n_rows, n_cols, -1);
    }

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

        free(matrix);
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

    if (error != MPI_SUCCESS)
        printf("Error %d\n", error);
    else
        printf("Success!\n");

    return 0;
}


void multiply_rowwise(double* local_matr, double* vector, long int n_rows, long int n_cols, double* result){
    for (long int i = 0; i < n_rows; ++i){
        double sum = 0;
        for (long int j = 0; j < n_cols; ++j){
            sum += local_matr[i * n_cols + j] * vector[j]; 
        }
        result[i] = sum;
    }

    return;
}



int main(int argc, char** argv){
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


    }

    long int local_n = n_rows / comm_sz;

    double* local_matr = malloc(local_n * n_cols * sizeof(double));
    double* vector = malloc(n_cols * sizeof(double));
    double* local_result = (double*) malloc(local_n * sizeof(double));
    double* result = (double*) malloc(n_rows * sizeof(double));

    int error;

    error = read_matrix_from_file(local_matr, local_n, n_rows, n_cols, my_rank, comm_sz, MPI_COMM_WORLD);
    if (error == -1){
        printf("Unable to locate matrix file '%s'\n", build_matrix_filename(n_rows, n_cols));
        return 0;
    }
    // print_matr(local_matr, local_n, n_cols, my_rank);

    error = read_vector_from_file(n_cols, my_rank, comm_sz, MPI_COMM_WORLD, vector);
    if (error == -1){
        printf("Unable to locate vector file '%s'\n", build_vector_filename(n_cols));
        return 0;
    }
    // print_vec(vector, n_cols, my_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    multiply_rowwise(local_matr, vector, local_n, n_cols, local_result);
    // print_vec(local_result, local_n, my_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    double local_elapsed = finish - start;
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, MAIN_PROCESS, MPI_COMM_WORLD);

    MPI_Gather(local_result, local_n, MPI_DOUBLE, result, local_n, MPI_DOUBLE, MAIN_PROCESS, MPI_COMM_WORLD);

    if (my_rank == MAIN_PROCESS){
        print_vec(result, n_rows, -1);
        printf("Elapsed time: %f s\n", elapsed);
    }

    free(local_matr);
    free(result);
    free(local_result);
    free(vector);

    MPI_Finalize();

    return 0;
}

