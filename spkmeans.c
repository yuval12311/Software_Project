#include <math.h> 
#include <float.h>
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
            printf("%.4f", a[i][j]);
            printf(j < n -1 ? ", " : "\n");
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
    off = -a[i][i]*a[i][i] - a[j][j]*a[j][j]; /* the difference in off^2 is just these two elements, since the sum of all elements squared is perserved under the rotation */
    a[i][j] = a[j][i] = 0;
    return off;
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

double** wam(vector* vectors, int n) {
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
        for (j = 0; j < n; ++j) {
            l[i][j] = - l[i][j];
        }
        l[i][i] += 1;
    }
    
    return l;
}


double** jacobi(vector* vectors, int n) {
    double **v, **p, **a, **temp, **ret;
    int imax, jmax, iterations = 0, l;
    double eps = 1.0e-15, c, s, t, off;
    double* eigenvals;


    a = lnorm(vectors, n);
    v = build_rot(0, 1, 0, 1, n);
    do {
        l = largest_off_diag(a, n);
        jmax = l%n; imax = l/n;
        t = compute_t(imax, jmax, a);
        c = compute_c(t);
        s = t * c;
        p = build_rot(imax, jmax, s, c, n);

        temp = v;
        v = mult(v, p, n);
        free_matrix(p);
        free_matrix(temp);
        off = update_a(imax, jmax, s, c, a, n);
        iterations++;
    } while(fabs(off) > eps && iterations <= 100 );

    eigenvals = diag(a, n);
    ret = malloc((n+1)*sizeof(double*));
    ret[0] = eigenvals;
    for (l = 0; l<n ;++l) {
        ret[1 + l] = v[l];
    }
    free(v);
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

int main() {

    return 0;


}

