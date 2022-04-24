#include <math.h> 
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"
#ifndef spkmeans
#define spkmeans
typedef struct {
    double** vecs;
    double* vals;
} eigen;

void free_matrix(double** a);
double** new_matrix(int n);
void freeVectors(vector* vectors, int k);
void freeVector(vector* v);
double** wam_c(vector* vectors, int n);
double** ddg_c(vector* vectors, int n);
double** lnorm_c(vector* vectors, int n);
eigen jacobi_c(double** a, int n);
#endif
