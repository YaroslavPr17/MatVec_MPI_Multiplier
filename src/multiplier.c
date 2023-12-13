#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LENGTH 128
#define MAIN_PROCESS 0


char* build_filename(int n_rows, int n_cols){
    char* filename = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    sprintf(filename, "%d_%d.txt", n_rows, n_cols);

    return filename;
}

void print_matr(double* matr, int n_rows, int n_cols, int proc_num){
    printf("%d\tMATRIX:\n", proc_num);
    for (int i = 0; i < n_rows; ++i){
        for (int j = 0; j < n_cols; ++j){
            printf("%lf ", matr[i * n_cols + j]);
        }   
        printf("\n");
    }

}


int read_matrix_from_file(double* local_matr[], int local_n, int n, int n_rows, int n_cols, int my_rank, int comm_sz, MPI_Comm comm){

    char* body = build_filename(n_rows, n_cols);
    char filename[MAX_FILENAME_LENGTH] = "./data/";
    strcat(filename, body);
    free(body);
    printf("Reading matrix from file '%s'\n", filename);

    double* matrix = (double*) malloc(n_rows * n_cols * sizeof(double));

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

    // int error;

    // MPI_Aint* displacements = (MPI_Aint*) malloc(sizeof(int));
    // printf("displacements:\n");
    // for (int i = 0; i < comm_sz; ++i){
    //     displacements[i] = &matrix[i][0] - &matrix[0][0];
    //     printf("%ld ", displacements[i]);
    // }
    // printf("\n");

    // int* blocklengths = (int*) malloc(sizeof(int));
    // printf("blocklengths:\n");
    // for (int i = 0; i < comm_sz; ++i){
    //     blocklengths[i] = n_cols;
    //     printf("%d ", blocklengths[i]);
    // }
    // printf("\n");


    // MPI_Datatype MPI_VECTOR;
    // error = MPI_Type_hindexed(
    //     n_rows,
    //     blocklengths,
    //     displacements,
    //     MPI_DOUBLE,
    //     &MPI_VECTOR
    // );
    // if (error != MPI_SUCCESS)
    //     printf("Error %d\n", error);
    // else
    //     printf("Success!\n");

    // error = MPI_Type_commit(&MPI_VECTOR);
    // if (error != MPI_SUCCESS)
    //     printf("Error %d\n", error);
    // else
    //     printf("Success!\n");

    // printf("Comm_size = %d\n", comm_sz);

    // error = MPI_Scatter(
    //     &matrix[0][0],
    //     local_n,
    //     MPI_VECTOR,
    //     local_matr,
    //     local_n,
    //     MPI_VECTOR,
    //     MAIN_PROCESS,
    //     comm
    // );
    // if (error != MPI_SUCCESS)
    //     printf("Error %d\n", error);
    // else
    //     printf("Success!\n");

    // MPI_Type_free(&MPI_VECTOR);

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


    read_matrix_from_file(local_matr, local_n, 4, n_rows, n_cols, my_rank, comm_sz, MPI_COMM_WORLD);

    // if (my_rank)
    //     print_matr(local_matr, 4, 3, my_rank);

    MPI_Finalize();

    return 0;
}

