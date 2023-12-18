# !/bin/bash

TYPE=$1

for n_proc in 1 2 6 12 24
do
    echo $n_proc
    for n_rows in 600 1800 3000 4200 5400 6600 7800 9000 10200
    do
        mpicc ./src/multiplier_$TYPE.c ./src/matr_utils.c ./src/utils.c -o ./out/multiplier -Wall -lm
        mpiexec -n $n_proc ./out/multiplier $n_rows $n_rows
    done
done
