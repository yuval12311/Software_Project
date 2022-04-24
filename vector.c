#include "vector.h"
void error() {
    printf("An Error Has Occurred\n");
    exit(1);
}

void input_error() {
    printf("Invalid Input!\n");
    exit(1);
}

void freeVector(vector* v) {
    free(v->data);
    free(v);
}

void freeVectors(vector* vectors, int k) {
    int i;
    for (i=0; i<k; ++i) {
        free(vectors[i].data);
    }
    free(vectors);
}

void printVector(vector* v) {
    int i;
    for (i=0; i<v->size; ++i) {
        printf("%f", v->data[i]);
        printf(i<v->size -1 ? ",\t" : "\n");
    }
}

void printVectors(vector* vs, int k) {
    int i;
    for (i=0; i<k; ++i) {
        printVector(vs+i);
    }
}
/* squared euclidean distance between vectors */
double sqr_dist(vector* v1, vector* v2) {
    double dist = 0;
    int i = 0;
    for (i = 0; i < v1->size; ++i) {
        dist += (v1->data[i] - v2->data[i])*(v1->data[i] - v2->data[i]);
    }
    return dist;
}

/* Adds vector b to vector a, in place */
void add(vector* a, vector* b) {   
    int i = 0; 
    for (i = 0; i < a->size; ++i) {
        a->data[i] = a->data[i] + b->data[i];
    }
}




vector* mult(double a, vector* v) {
    
    vector* v2 = (vector*)malloc(sizeof(vector));
    int i = 0;
    if (v2 == NULL) error();
    v2->size = v->size;
    
    v2->data = (double*)malloc(sizeof(double) * v2->size);
    if (v2->data == NULL) error();

    for (i = 0; i < v2->size; ++i) {
        v2->data[i] = a* v->data[i];
    }
    return v2;
}