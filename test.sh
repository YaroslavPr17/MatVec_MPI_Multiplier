# !/bin/bash

TYPE=$1

N_ROWS_START_VALUE=120
N_ROWS_FINISH_VALUE=1201
N_ROWS_STEP=120

N_COLS=60000

for n_proc in 1 2 6 12 24
do
    echo $n_proc
    for ((n_rows=$N_ROWS_START_VALUE; n_rows < $N_ROWS_FINISH_VALUE; n_rows += $N_ROWS_STEP))
    do
        mpicc ./src/multiplier_$TYPE.c ./src/matr_utils.c ./src/utils.c -o ./out/multiplier -Wall -lm
        mpiexec -n $n_proc ./out/multiplier $n_rows $N_COLS
    done
done
