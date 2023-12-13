#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>


#define MAX_SHORT_NAME_LENGTH 64
#define MAX_FILENAME_LENGTH 128

double generate_random_sign(){
    return pow(-1, (double)(rand() % 10));
}


int main(int argc, char** argv){    
    const long MAX_VAL = LONG_MAX;

    srand(time(NULL));

    short is_vector = 0, is_matrix = 0;

    char* type = argv[1];
    if (!strcmp(type, "vector"))
        is_vector = 1;
    else if (!strcmp(type, "matrix"))
        is_matrix = 1;
    else{
        printf("Unknown data type. Exitting...");
    }
    

    long n_rows = strtol(argv[2], NULL, 10);
    long n_cols = -1;
    if (is_matrix)
        n_cols = strtol(argv[3], NULL, 10);

    char filename[MAX_FILENAME_LENGTH];
    char short_filename[MAX_SHORT_NAME_LENGTH];
    
    if (is_matrix){
        sprintf(short_filename, "%s_%ld_%ld.txt", type, n_rows, n_cols);
    }
    else {
        sprintf(short_filename, "%s_%ld.txt", type, n_rows);
    }

    strcpy(filename, "./data/");
    strcat(filename, short_filename);

    printf("---> %s will be generated with size (%ld, %ld).\n", type, n_rows, n_cols);

    FILE *fp = fopen(filename, "w");

    if (fp == NULL){
        printf("File opening error! (%s)\n", filename); 
        return 0;
    }

    if (is_matrix)
        for (long i = 0; i < n_rows; ++i){
            for (long j = 0; j < n_cols; ++j){
                fprintf(fp, "%.3lf ", rand() % MAX_VAL * generate_random_sign() / 10e3);
            }
            fprintf(fp, "\n");
        }
    else
        for (long i = 0; i < n_rows; ++i){
            fprintf(fp, "%.3lf\n", rand() % MAX_VAL * generate_random_sign() / 10e3);
        }

    fclose(fp);

    return 0;
}
