#ifndef MATR_UTILS_H
#define MATR_UTILS_H

void build_matrix_filename(long int n_rows, long int n_cols, char* filename);
void build_vector_filename(long int n_elems, char* filename);
void print_matr(double* matr, long int n_rows, long int n_cols, int proc_num);
void print_vec(double* vec, long int n_elems, int proc_num);
int load_matr(long int n_rows, long int n_cols, double* matrix);
int load_vec(long int n_rows,  double* vector);


#endif /* MATR_UTILS_H */
