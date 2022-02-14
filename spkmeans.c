#include <math.h> 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int size;
    double* data;
} vector;

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

void print_matrix(double** a, int n) {
    int i, j;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            pritnf("%.4f", a[i][j]);
            printf(i < n -1 ? ",\t" : "\n");
        }
    }
    
}

/* returns new square matrix, initialized with 0's */
double** new_matrix(int n) {
    double** a;
    double* p;
    int i;

    a = malloc (sizeof(double*) * n);
    if (a == NULL) error();

    p = calloc (n*n, sizeof(double));
    if (a == NULL) error();

    for (i = 0; i < n; ++i) {
        a[i] = p + i*n;
    }

    return a;
}

/* frees a matrix */
void free_matrix(double** a) {
    free(a[0]);
    free(a);
}

/* Multiplies to n*n matrices */
double** mult(double** a, double** b, int n) {
    double** c;
    int i, j, k;

    c = new_matrix(n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            for (k = 0; k < n; ++k){
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    return c;
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

/* inverse square root of a diagonal matrix, inplace */
void inv_sqrt(double** d, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        d[i][i] = 1 / sqrt(d[i][i]);
    }
}


double** wam(vector* vectors, int n) {
    double** w;
    int i,j;

    w = new_matrix(n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            w[i][j] = exp(-sqrt(sqr_dist(vectors + i, vectors + j))/2);
        }
    }
    return w;
}

double** ddg_from_w(double** w, int n) {
    double** d;
    int i;
    int j;

    d = new_matrix(n);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            d[i][i] += w[i][j];
        }
    }
    return d;
}

double** ddg(vector* vectors, int n) {
    double **w, **d;
    w = wam(vectors, n);
    d = ddg_from_w(w, n);
    free_matrix(w);
    return d;
}

double** lnorm(vector* vectors, int n) {
    double **w, **d, **l , **temp;
    int i,j;

    w = wam(vectors, n);
    d = ddg_from_w(w, n);
    inv_sqrt(d, n);

    l = mult(d, w, n);
    temp = l;
    l = mult(l, d, n);
    free_matrix(temp);
    free_matrix(d);
    free_matrix(w);

    for (i = 0; i < n; ++i) {
        l[i][i] = 1 - l[i][i];
    }
    
    return l;
}


