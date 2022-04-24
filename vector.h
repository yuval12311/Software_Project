#include <math.h> 
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef vec
#define vec
typedef struct vector {
    int size;
    double* data;
} vector;
#endif
void error();
void input_error();
void freeVector(vector* v);
void freeVectors(vector* vectors, int k);
void printVector(vector* v);
void printVectors(vector* vs, int k);
double sqr_dist(vector* v1, vector* v2);
void add(vector* a, vector* b);
vector* mult(double a, vector* v);
