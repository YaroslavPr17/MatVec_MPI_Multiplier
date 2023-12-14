#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LENGTH 128
#define MAIN_PROCESS 0

void process_error(int err){
    if (err != MPI_SUCCESS){
        printf("Error %d\n", err);

        char* str = malloc(128 * sizeof(char));
        int str_len; 

        MPI_Error_string(err, str, &str_len);

        printf("Error: %s, length: %d\n", str, str_len);
    }
}

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


int read_nums_from_file(long int n_elems, int local_n, int my_rank, int comm_sz, MPI_Comm comm, double* nums){

    double* vector;

    if (my_rank == MAIN_PROCESS){
        char* body = build_vector_filename(n_elems);
        char filename[MAX_FILENAME_LENGTH] = "./data/";
        strcat(filename, body);
        // printf("Before free\n");
        free(body);
        printf("Reading vector from file '%s'...\n", filename);

        vector = (double*) malloc(n_elems * sizeof(double));

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


    error = MPI_Scatter(
        vector,
        local_n,
        MPI_DOUBLE,
        nums,
        local_n,
        MPI_DOUBLE,
        MAIN_PROCESS,
        comm
    );

    process_error(error);

    // if (my_rank == MAIN_PROCESS)
    //     free(vector);

    return 0;
}


void print_matr(double* matr, long int n_rows, long int n_cols, int proc_num){
    printf("%d\tMATRIX (%ld x %ld):\n", proc_num, n_rows, n_cols);
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

    // printf("local_n = %ld\n", local_n);
    // printf("n_rows = %ld\n", n_rows);
    // printf("n_cols = %ld\n", n_cols);
    // printf("my_rank = %d\n", my_rank);
    // printf("comm_sz = %d\n", comm_sz);


    if (my_rank == MAIN_PROCESS){
        char* body = build_matrix_filename(n_rows, n_cols);
        char filename[MAX_FILENAME_LENGTH] = "./data/";
        strcat(filename, body);
        // printf("Before free\n");
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
                15,
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
            15, 
            MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE
        );
        process_error(error);
        // print_matr(local_matr, n_rows, local_n, my_rank);
    }

    printf("AFTER TRANSPORTATION\n");

    // free(matrix);
    

    MPI_Type_free(&mpi_vec);

    process_error(error);

    return 0;
}


void multiply_colwise(double* local_matr, double* nums, long int n_rows, long int local_n, int my_rank, int comm_sz, double* result){
    
    for (long int j = 0; j < local_n; ++j){
        for (long int i = 0; i < n_rows; ++i){
            local_matr[i * local_n + j] *= nums[j]; 
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


    // double* matrix = (double*) malloc(n_rows * n_cols * sizeof(double));

    // MPI_Gather(local_matr, n_rows * local_n, MPI_DOUBLE, matrix, n_rows * local_n, MPI_DOUBLE, MAIN_PROCESS, MPI_COMM_WORLD);

    // if (my_rank == MAIN_PROCESS)
    //     print_matr(matrix, n_rows, n_cols, -1);

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
        if (n_cols % comm_sz != 0){
            printf("\nERROR!!!\n%ld mod %d = %ld. Unable to parallellize task.\n", n_rows, comm_sz, n_rows % comm_sz);
            return 0;
        }
    }

    long int local_n = n_cols / comm_sz;

    double* local_matr = malloc(local_n * n_rows * sizeof(double));
    double* nums = malloc(local_n * sizeof(double));
    double* local_result = (double*) malloc(local_n * sizeof(double));
    double* result = (double*) malloc(n_rows * sizeof(double));

    int error;

    error = read_matrix_from_file(local_matr, local_n, n_rows, n_cols, my_rank, comm_sz, MPI_COMM_WORLD);
    if (error == -1){
        printf("Unable to locate matrix file '%s'\n", build_matrix_filename(n_rows, n_cols));
        return 0;
    }
    // if (my_rank == comm_sz - 1)
    //     print_matr(local_matr, n_rows, local_n, my_rank);

    error = read_nums_from_file(n_cols, local_n, my_rank, comm_sz, MPI_COMM_WORLD, nums);
    if (error == -1){
        printf("Unable to locate initial vector file '%s'\n", build_vector_filename(n_cols));
        return 0;
    }
    // if (my_rank == comm_sz - 1)
    //     print_vec(nums, local_n, my_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    printf("BEFORE MULTIPLICATON\n");

    multiply_colwise(local_matr, nums, n_rows, local_n, my_rank, comm_sz, result);

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
    free(result);
    free(local_result);
    free(nums);

    MPI_Finalize();

    return 0;
}

