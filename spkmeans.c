#include <math.h> 
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"
typedef struct {
    double** vecs;
    double* vals;
} eigen;




void print_arr_jacobi(double* a, int n) {
    int j;
    for (j = 0; j < n; ++j) {
            if (fabs(a[j]) >= 0.00005) printf("%.4f", a[j]);
            else printf("0.0000");
            printf(j < n -1 ? "," : "\n");
    }
}
void print_arr(double* a, int n) {
    int j;
    for (j = 0; j < n; ++j) {
            printf("%.4f", a[j]);
            
            printf(j < n -1 ? "," : "\n");
    }
}

void print_matrix(double** a, int n) {
    int i;

    for (i = 0; i < n; ++i) {
        print_arr(a[i], n);
    }
    
}

void print_matrix_jacobi(double** a, int n) {
    int i;

    for (i = 0; i < n; ++i) {
        print_arr_jacobi(a[i], n);
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

double** vectors_to_matrix(vector* vectors, int n) {
    double** matrix;
    int i,j;    
    matrix = new_matrix(n);
    for (i = 0; i < n; ++i) {
        for (j = i; j < n; ++j) {
            if (vectors[i].data[j] != vectors[j].data[i]) input_error();
            matrix[i][j] = matrix[j][i] = vectors[i].data[j];
        }
    }
    return matrix;
}

/* Multiplies to n*n matrices */
double** mat_mult(double** a, double** b, int n) {
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

int in(char* str, char* strs[], int size) {
    int i;
    for (i = 0; i < size; ++i) {
        if (strcmp(str, strs[i]) == 0) return i;
    }
    return size;
}


/* inverse square root of a diagonal matrix, inplace */
void inv_sqrt(double** d, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        d[i][i] = 1.0 / sqrt(d[i][i]);
    }
}


double offcalc(double** a, int n) {
    int i, j;
    double sum = 0;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i != j) sum += a[i][j]*a[i][j];
        }
    }
    return sum;
}


int largest_off_diag(double** a, int n) {
    int i,j,imax=0,jmax=1;
    for (i = 0; i < n; ++i) {
        for( j = i + 1; j < n; ++j) {
            if (fabs(a[i][j]) > fabs(a[imax][jmax])) {
                imax = i;
                jmax = j;
            }
        }
    }

    return imax * n + jmax;

}

int sign(double x) {
    return x >= 0 ? 1 : -1;
}

double compute_t(int i, int j, double** a) {
    double theta = (a[j][j] - a[i][i])/(2*a[i][j]);
    return sign(theta) / (fabs(theta) + sqrt(theta * theta + 1));
}
double compute_c(double t) {
    return 1/sqrt(t*t+1);
}


double** build_rot(int i, int j, double s, double c, int n) {
    double** p;
    int k;
    
    
    p = new_matrix(n);
    for (k = 0; k < n; ++k) {p[k][k] = 1;}
    p[i][i] = c;
    p[j][j] = c;
    p[i][j] = s;
    p[j][i] = -s;

    return p;
}


/* updates A according to the rotation matrix P, returns the difference in off^2 between the old and the new matrix */
double update_a(int i, int j, double s, double c, double** a, int n) {
    int r;
    double off, temp1, temp2;
    for(r = 0; r < n; ++r) {
        if (r != i && r != j) {
            temp1 = a[r][i];
            temp2 = a[r][j];
            a[r][i] = a[i][r] = c*temp1 - s*temp2;
            a[r][j] = a[j][r] = c*temp2 + s*temp1;
        }
    }
    temp1 = a[i][i]; temp2 = a[j][j];
    a[i][i] = c*c*temp1 + s*s*temp2 - 2*s*c*a[i][j];
    a[j][j] = s*s*temp1 + c*c*temp2 + 2*s*c*a[i][j];
    off =temp1*temp1 -a[i][i]*a[i][i] + temp2*temp2 - a[j][j]*a[j][j]; /* the difference in off^2 is just these two elements, since the sum of all elements squared is perserved under the rotation */
    a[i][j] = a[j][i] = 0;
    return  off;
}

void update_v(int i, int j, double s, double c, double** v, int n) {
    int r;
    double temp1, temp2;
    for(r = 0; r < n; ++r) {
        temp1 = v[r][i]; temp2 = v[r][j];
        v[r][i] = temp1 * c - temp2 * s;
        v[r][j] = temp1 * s + temp2 * c;
    }
}

double* diag(double** a, int n) {
    double* diag;
    int i;
    diag = malloc(n*sizeof(double));
    if (diag == NULL) error();
    for (i = 0; i < n; ++i) {
        diag[i] = a[i][i];
    }
    return diag;
}

double** wam_c(vector* vectors, int n) {
    double** w;
    int i,j;

    w = new_matrix(n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            w[i][j] = i == j ? 0 : exp(-sqrt(sqr_dist(vectors + i, vectors + j))/2);
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

double** ddg_c(vector* vectors, int n) {
    double **w, **d;
    w = wam_c(vectors, n);
    d = ddg_from_w(w, n);
    free_matrix(w);
    return d;
}

double** lnorm_c(vector* vectors, int n) {
    double **w, **d, **l , **temp;
    int i,j;

    w = wam_c(vectors, n);
    d = ddg_from_w(w, n);
    inv_sqrt(d, n);

    l = mat_mult(d, w, n);
    temp = l;
    l = mat_mult(l, d, n);
    free_matrix(temp);
    free_matrix(d);
    free_matrix(w);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            l[i][j] = - l[i][j];
        }
        l[i][i] += 1;
    }
    
    return l;
}


eigen jacobi_c(double** a, int n) {
    double **v;
    int imax, jmax, iterations = 1, l;
    double eps = 1.0e-5, c, s, t, off=0;
    double* eigenvals;
    eigen ret;

    
    v = build_rot(0, 1, 0, 1, n);
     do {
        l = largest_off_diag(a, n);
        jmax = l%n; imax = l/n;
        if (a[imax][jmax] == 0) break;
        t = compute_t(imax, jmax, a);
        c = compute_c(t);
        s = t * c;
        /* p = build_rot(imax, jmax, s, c, n);

        temp = v;
        v = mat_mult(v, p, n);
        free_matrix(p);
        free_matrix(temp); */
        update_v(imax, jmax, s, c, v, n);
        off = update_a(imax, jmax, s, c, a, n);
        /* print_matrix(v, n);
        printf("\n");
        print_matrix(a, n);
        printf("\n%d, %17g\n\n", iterations, -off); */
        iterations++;
    }while(-off > eps && iterations <= 100 );

    eigenvals = diag(a, n);
    
    free_matrix(a);

    ret.vecs = v;
    ret.vals = eigenvals;
    return ret;
}

int countLines(char* fname) {
    FILE *ifp = fopen(fname, "r");
    int i = 0;
    char c;
    if (ifp == NULL) error();

    while ((c=fgetc(ifp)) != EOF) {
        if (c == '\n') ++i;
    }
    fclose(ifp);
    return i;
}

int vectorSize(char* fname) {
    FILE *ifp = fopen(fname, "r");
    int i = 1;
    char c;
    if (ifp == NULL) error();
    
    while ((c=fgetc(ifp)) != '\n') {
        if (c == ',') ++i;
    }
    fclose(ifp);
    return i;
}

vector* read(char* fname, int n) {
    int d;
    int i;
    int j;
    vector* vectors;
    FILE *ifp;
    vectors = (vector*) calloc(n, sizeof(vector));
    if (vectors == NULL) error();
    d = vectorSize(fname);
    i = 0;
    ifp = fopen(fname, "r");
    if (ifp == NULL) error();
    while(i < n) {
        (vectors+i)->size = d;
        (vectors+i)->data = (double*) calloc(d, sizeof(double));
        if ((vectors+i)->data == NULL) error();
        for (j = 0; j < d; ++j) {
            fscanf(ifp, "%lf", (vectors+i)->data + j);
            fgetc(ifp);
        }
        ++i;
    }
    fclose(ifp);
    return vectors;
}

int main(int argc, char* argv[]) {
    char* file;
    double** matrix,** matrix2;
    int type, n;
    vector* vectors;
    char* functions[] = {"wam", "ddg", "lnorm", "jacobi"};
    eigen ret;
    
    if (argc != 3) input_error();
    file = argv[2];
    type = in(argv[1], functions, 4);
    if (type > 3) input_error();
    n = countLines(file);
    vectors = read(file, n);
    switch (type)
    {
    case 0:
        matrix = wam_c(vectors, n);
        break;
    case 1:
        matrix = ddg_c(vectors, n);
        break;
    case 2: 
        matrix = lnorm_c(vectors, n);
        break;
    case 3:
        if (vectors->size != n) input_error();
        matrix2 = vectors_to_matrix(vectors, n);
        ret = jacobi_c(matrix2, n);
        matrix = ret.vecs;
        
        print_arr(ret.vals, n);
        free(ret.vals);
        
        break;
    default:
        matrix = new_matrix(1);
        break;
    }
    if (type!=3) print_matrix(matrix, n);
    else print_matrix_jacobi(matrix, n);
    
    free_matrix(matrix);
    freeVectors(vectors, n);
    
    /* double** p = build_rot(0, 1, 0, 1, 5);
    free_matrix(p); */
    return 0;


}

