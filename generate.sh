# !/bin/bash

TYPE=$1

N_ROWS_START_VALUE=$2
N_ROWS_FINISH_VALUE=$3
N_ROWS_STEP=$4

N_COLS_START_VALUE=$5
N_COLS_FINISH_VALUE=$6
N_COLS_STEP=$7

if [[ $TYPE = "matrix" ]]; then 
    for ((n_rows=$N_ROWS_START_VALUE; n_rows < $N_ROWS_FINISH_VALUE; n_rows += $N_ROWS_STEP))
    do
        for ((n_cols=$N_COLS_START_VALUE; n_cols < $N_COLS_FINISH_VALUE; n_cols += $N_COLS_STEP))
        do
            ./out/generator $TYPE $n_rows $n_cols
            sleep 1
        done
    done
else
    if [[ $TYPE = "vector" ]]; then
        for ((n_rows=$N_ROWS_START_VALUE; n_rows < $N_ROWS_FINISH_VALUE; n_rows += $N_ROWS_STEP))
        do
            ./out/generator $TYPE $n_rows
            sleep 1
        done
    else
        echo "Unknown type of object to generate. Exitting..."
    fi
fi

