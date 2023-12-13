#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LENGTH 128
#define MAIN_PROCESS 0


char* build_matrix_filename(int n_rows, int n_cols){
    char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    sprintf(filename, "matrix_%d_%d.txt", n_rows, n_cols);

    return filename;
}

char* build_vector_filename(int n_elems){
    char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    sprintf(filename, "vector_%d.txt", n_elems);

    return filename;
}


void print_vec(double* vec, int n_elems, int proc_num){
    printf("%d\tVECTOR:\n", proc_num);
    for (int i = 0; i < n_elems; ++i){
        printf("%lf\n", vec[i]);
    }
    return;
}


int read_vector_from_file(int n_elems, int my_rank, int comm_sz, MPI_Comm comm, double* vector){

    if (my_rank == MAIN_PROCESS){
        char* body = build_vector_filename(n_elems);
        char filename[MAX_FILENAME_LENGTH] = "./data/";
        strcat(filename, body);
        free(body);
        printf("Reading vector from file '%s'\n", filename);

        FILE *fp = fopen(filename, "r");
        if (fp == NULL){ 
            return -1;
        }

        for (int i = 0; i < n_elems; ++i){
            fscanf(fp, "%lf", &vector[i]); 
        }

        print_vec(vector, n_elems, -1);
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


void print_matr(double* matr, int n_rows, int n_cols, int proc_num){
    printf("%d\tMATRIX:\n", proc_num);
    for (int i = 0; i < n_rows; ++i){
        for (int j = 0; j < n_cols; ++j){
            printf("%lf ", matr[i * n_cols + j]);
        }   
        printf("\n");
    }
    return;
}


int read_matrix_from_file(double* local_matr, int local_n, int n_rows, int n_cols, int my_rank, int comm_sz, MPI_Comm comm){
    double* matrix;

    if (my_rank == MAIN_PROCESS){
        char* body = build_matrix_filename(n_rows, n_cols);
        char filename[MAX_FILENAME_LENGTH] = "./data/";
        strcat(filename, body);
        free(body);
        printf("Reading matrix from file '%s'\n", filename);

        matrix = (double*) malloc(n_rows * n_cols * sizeof(double));

        FILE *fp = fopen(filename, "r");
        if (fp == NULL){ 
            return -1;
        }

        for (int i = 0; i < n_rows; ++i){
            for (int j = 0; j < n_cols; ++j){
                fscanf(fp, "%lf", &matrix[i * n_cols + j]);
            }   
        }

        print_matr(matrix, n_rows, n_cols, -1);
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



int main(int argc, char** argv){
    int n_rows = strtol(argv[1], NULL, 10);
    int n_cols = strtol(argv[2], NULL, 10);

    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int local_n = n_rows / comm_sz;

    double* local_matr = malloc(local_n * n_cols * sizeof(double));
    double* vector = malloc(n_cols * sizeof(double));

    read_matrix_from_file(local_matr, local_n, n_rows, n_cols, my_rank, comm_sz, MPI_COMM_WORLD);
    print_matr(local_matr, local_n, n_cols, my_rank);

    read_vector_from_file(n_cols, my_rank, comm_sz, MPI_COMM_WORLD, vector);
    // print_vec(vector, n_cols, my_rank);




    free(local_matr);
    free(vector);

    MPI_Finalize();

    return 0;
}

