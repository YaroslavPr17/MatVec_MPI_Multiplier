#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matr_utils.h"
#include "utils.h"
#include "constants.h"


void distribute_data (double* matrix, double* vector, long int n_rows, long int n_cols, long int local_n, int my_rank, int comm_sz, MPI_Comm comm, double* local_matr, double* local_vec){

    MPI_Datatype mpi_vec;

    MPI_Type_vector(
        n_rows,
        local_n,
        n_cols,
        MPI_DOUBLE,
        &mpi_vec
    );

    MPI_Type_commit(&mpi_vec);

    int error;

    // double* packed_matrix = (double*) malloc(n_rows * n_cols * sizeof(double));

    // int position = 0;

    // int pack_size;
    // MPI_Pack_size(local_n, mpi_vec, MPI_COMM_WORLD, &pack_size);


    // if (my_rank == MAIN_PROCESS){
    //     // printf("Packing...\n");
    //     // MPI_Pack(
    //     //     &matrix[0],
    //     //     local_n,
    //     //     mpi_vec,
    //     //     packed_matrix,
    //     //     pack_size,
    //     //     &position,
    //     //     MPI_COMM_WORLD
    //     // );

    //     printf("Sending...\n");
    //     error = MPI_Scatter(
    //         &matrix[0],
    //         1,
    //         mpi_vec,
    //         local_matr,
    //         1,
    //         mpi_vec,
    //         MAIN_PROCESS,
    //         comm
    //     );

    // }
    // else {
    //     void* recv_buffer = malloc(n_rows * n_cols * sizeof(double));

    //     printf("Receiving...\n");
    //     error = MPI_Scatter(
    //         NULL,
    //         -1,
    //         MPI_PACKED,
    //         &recv_buffer,
    //         1,
    //         MPI_PACKED,
    //         MAIN_PROCESS,
    //         comm
    //     );
        
    //     printf("Unpacking...\n");
    //     MPI_Unpack(
    //         recv_buffer,
    //         pack_size,
    //         &position,
    //         &local_matr[0],
    //         local_n,
    //         mpi_vec,
    //         MPI_COMM_WORLD
    //     );

    //     print_matr(local_matr, n_rows, local_n, my_rank);

    // }

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("All here!\n");
    
    printf("BEFORE TRANSPORTATION\n");

    if (my_rank == MAIN_PROCESS){

        // double* packed_matrix = (double*) malloc(n_rows * local_n * sizeof(double));

        int position;

        int pack_size;
        MPI_Pack_size(1, mpi_vec, MPI_COMM_WORLD, &pack_size);

        for (int i = 1; i < comm_sz; ++i){

            position = 0;

            // printf("Packing (%d)...\n", i);
            error = MPI_Pack(
                &matrix[i * local_n],
                1,
                mpi_vec,
                local_matr,
                pack_size,
                &position,
                MPI_COMM_WORLD
            );
            process_error(error);

            // printf("Sending (%d)...\n", i);
            error = MPI_Send(
                local_matr,
                n_rows * local_n,
                MPI_DOUBLE,
                i,
                SUBMATR_TAG,
                MPI_COMM_WORLD
            );
            process_error(error);
        }

        // print_matr(matrix, n_rows, n_cols, -1);

        position = 0;

        error = MPI_Pack(
            &matrix[0],
            1,
            mpi_vec,
            local_matr,
            pack_size,
            &position,
            MPI_COMM_WORLD
        );
        process_error(error);

        // print_matr(local_matr, n_rows, local_n, -100);

    }
    else {
        // printf("Receiving (%d)...\n", my_rank);
        error = MPI_Recv( 
            &local_matr[0],
            n_rows * local_n, 
            MPI_DOUBLE,
            MAIN_PROCESS,
            SUBMATR_TAG, 
            MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE
        );
        process_error(error);
        // print_matr(local_matr, n_rows, local_n, my_rank);
    }


    error = MPI_Scatter(
        vector,
        local_n,
        MPI_DOUBLE,
        local_vec,
        local_n,
        MPI_DOUBLE,
        MAIN_PROCESS,
        comm
    );
    process_error(error);


    printf("AFTER TRANSPORTATION\n");
    

    error = MPI_Type_free(&mpi_vec);
    process_error(error);

    return;
}


void multiply_colwise(double* local_matr, double* local_vec, long int n_rows, long int local_n, int my_rank, int comm_sz, double* result){
    
    for (long int j = 0; j < local_n; ++j){
        for (long int i = 0; i < n_rows; ++i){
            local_matr[i * local_n + j] *= local_vec[j]; 
        }
    }

    // if (my_rank == comm_sz - 1){
    //     print_matr(local_matr, n_rows, local_n, my_rank);
    //     // printf("n_rows = %ld, local_n = %ld\n", n_rows, local_n);
    // }

    long int n_cols = local_n * comm_sz;

    double* columns = (double*) malloc(n_cols * sizeof(double));
    for (long int i = 0; i < n_rows; ++i){
        double sum = 0.0;
        for (long int j = 0; j < local_n; ++j){
            sum += local_matr[i * local_n + j];
        }
        columns[i] = sum;
    } 

    // if (my_rank == comm_sz - 1)
    //     print_vec(columns, n_rows, my_rank);   

    MPI_Reduce(         // Должны вызывать все процессы
        columns,        // Буфер, который отсылаем
        result,         // Буфер ответа. Важен только тот, который идёт во 0-й поток
        n_rows,         // Количество чисел
        MPI_DOUBLE,     // 
        MPI_SUM,        // Как агрегируем
        MAIN_PROCESS,              // Куда отправляем результат
        MPI_COMM_WORLD
    );

    free(columns);


    // double* matrix = (double*) malloc(n_rows * n_cols * sizeof(double));

    // MPI_Gather(local_matr, n_rows * local_n, MPI_DOUBLE, matrix, n_rows * local_n, MPI_DOUBLE, MAIN_PROCESS, MPI_COMM_WORLD);

    // if (my_rank == MAIN_PROCESS)
    //     print_matr(matrix, n_rows, n_cols, -1);

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
        if (n_cols % comm_sz != 0){
            printf("\nERROR!!!\n%ld mod %d = %ld. Unable to parallellize task.\n", n_cols, comm_sz, n_rows % comm_sz);
            return 0;
        }
    }

    long int local_n = n_cols / comm_sz;


    double* matrix;
    double* vector;
    double* result;


    if (my_rank == MAIN_PROCESS){
        printf("n_rows = %ld\n", n_rows);
        printf("n_cols = %ld\n", n_cols);
        printf("local_n = %ld\n", local_n);
        printf("comm_sz = %d\n", comm_sz);
        printf("my_rank = %d\n", my_rank);

        matrix = (double*) malloc(n_rows * n_cols * sizeof(double));
        vector = (double*) malloc(n_cols * sizeof(double));
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


    double* local_matr = malloc(local_n * n_rows * sizeof(double));
    double* local_vec = malloc(local_n * sizeof(double));


    distribute_data(matrix, vector, n_rows, n_cols, local_n, my_rank, comm_sz, MPI_COMM_WORLD, local_matr, local_vec);
    // if (my_rank == comm_sz - 1){
    //     print_matr(local_matr, n_rows, local_n, my_rank);
    //     print_vec(local_vec, local_n, my_rank);
    // } 

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    printf("BEFORE MULTIPLICATON\n");

    multiply_colwise(local_matr, local_vec, n_rows, local_n, my_rank, comm_sz, result);

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    double local_elapsed = finish - start;
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, MAIN_PROCESS, MPI_COMM_WORLD);

    if (my_rank == MAIN_PROCESS){
        printf("RESULT:\n");
        print_vec(result, n_rows, -1);
        printf("Elapsed time: %f s\n", elapsed);
    }


    free(local_matr);
    free(local_vec);
    if (my_rank == MAIN_PROCESS){
        free(result);
        free(vector);
        free(matrix);
    }

    MPI_Finalize();

    return 0;
}

