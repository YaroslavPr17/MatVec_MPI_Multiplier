#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matr_utils.h"
#include "constants.h"


void build_matrix_filename(long int n_rows, long int n_cols, char* filename){
    sprintf(filename, "matrix_%ld_%ld.txt", n_rows, n_cols);
    return;
}

void build_vector_filename(long int n_elems, char* filename){
    sprintf(filename, "vector_%ld.txt", n_elems);
    return;
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

void print_vec(double* vec, long int n_elems, int proc_num){
    printf("%d\tVECTOR:\n", proc_num);
    for (long int i = 0; i < n_elems; ++i){
        printf("%d\t%lf\n", proc_num, vec[i]);
    }
    return;
}


int load_matr(long int n_rows, long int n_cols, double* matrix){
    char* body = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    build_matrix_filename(n_rows, n_cols, body);
    char filename[MAX_FILENAME_LENGTH] = "./data/";
    strcat(filename, body);
    free(body);
    printf("Reading matrix from file '%s'...\n", filename);

    FILE *fp = fopen(filename, "r");
    if (fp == NULL){
        return -1;
    }

    for (long int i = 0; i < n_rows; ++i){
        for (long int j = 0; j < n_cols; ++j){
            fscanf(fp, "%lf", &matrix[i * n_cols + j]);
        }   
    }

    return 0;
}


int load_vec(long int n_rows,  double* vector){
    char* body = (char*) malloc(MAX_FILENAME_LENGTH * sizeof(char));
    build_vector_filename(n_rows, body);
    char filename[MAX_FILENAME_LENGTH] = "./data/";
    strcat(filename, body);
    free(body);
    printf("Reading vector from file '%s'...\n", filename);

    FILE *fp = fopen(filename, "r");
    if (fp == NULL){ 
        return -1;
    }

    for (long int i = 0; i < n_rows; ++i){
        fscanf(fp, "%lf", &vector[i]); 
    }

    return 0;
}