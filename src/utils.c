#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "constants.h"


void process_error(int err){
    if (err != MPI_SUCCESS){
        printf("Error %d\n", err);

        char* str = malloc(STR_DEFAULT_LENGTH * sizeof(char));
        int str_len; 

        MPI_Error_string(err, str, &str_len);

        printf("Error: %s, length: %d\n", str, str_len);

        free(str);
    }
}


void get_2_most_closest_multipliers(long int number, int* dividers){
    int sroot = (int) sqrt((double) number);
    for (int cur_div = sroot; cur_div > 0; --cur_div){
        if (number % cur_div == 0){
            dividers[0] = cur_div;
            dividers[1] = number / cur_div;

            return;
        }
    }
    return;
}

