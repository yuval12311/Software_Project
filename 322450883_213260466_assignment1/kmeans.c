#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define epsilon 0.001
int numVectors;


typedef struct {
    int size;
    double* data;
} vector;


/* Adds vector b to vector a, in place */
void add(vector* a, vector* b) {   
    int i = 0; 
    for (i = 0; i < a->size; ++i) {
        a->data[i] = a->data[i] + b->data[i];
    }
}

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

/* return euclidean distace squared between v1 and v2 */
double sqr_dist(vector* v1, vector* v2) {
    double dist = 0;
    int i = 0;
    for (i = 0; i < v1->size; ++i) {
        dist += (v1->data[i] - v2->data[i])*(v1->data[i] - v2->data[i]);
    }
    return dist;
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

void initalize_cluster_sizes(int* sizes, int k) {
    int i = 0;
    for (i = 0; i < k; ++i) {
        sizes[i] = 0;
    }
}

void initalize_cluster_sums(vector* sums, int k) {
    int i = 0;
    int j=0;
    for (i = 0; i < k; ++i) {
        for(j =0;j < sums[i].size; j++) {
            sums[i].data[j]=0;
        }
    }
}

int closest_centroid(vector* v, vector* centroids, int k) {
    int min_i = 0;
    double min_dist = sqr_dist(v, centroids);
    int boba;
    for (boba = 0; boba < k; ++boba) {
        if (min_dist > sqr_dist(v, centroids + boba)) {
            min_i = boba;
            min_dist = sqr_dist(v, centroids + boba);
        }
    }
   
    return min_i;

}

void copydata(vector* v1, vector* v2) {
    int i;
    v1->size = v2->size;
    for (i = 0; i<v1->size; ++i) {
        v1->data[i] = v2->data[i];
    }
}
vector* kmeans(vector* vectors, int k, int maxiter) {
    vector* mu = (vector*) calloc(k, sizeof(vector));
    vector* newmui;
    int* cluster_sizes = (int*) calloc(k, sizeof(int));
    vector* cluster_sums = (vector*) calloc(k, sizeof(vector));
    int bob=0;
    int finished = 0;
    int iter = 0;
    int i=0;
    int j=0;

    if (mu == NULL || cluster_sizes == NULL || cluster_sums == NULL) 
        error();

    for(bob = 0; bob < k;bob++)
    {
        cluster_sums[bob].size = vectors->size;
        cluster_sums[bob].data = calloc(vectors->size, sizeof(double)); 
        if (cluster_sums[bob].data == NULL) error(); 
    }
    for(bob = 0; bob < k; ++bob) {
        (mu+bob)->data =(double*) malloc(sizeof(double)*vectors->size);
        if ((mu+bob)->data == NULL) error(); 
        copydata(mu+bob, vectors+bob);
    }
    while (!finished) {
        initalize_cluster_sizes(cluster_sizes,k);
        initalize_cluster_sums(cluster_sums,k);


        for(i=0;i<numVectors;i++){
            vector* x = vectors+i;


            j = closest_centroid(x, mu, k);
            

            cluster_sizes[j] += 1;
            add(cluster_sums+j, x);
        }
        
        iter++;
        finished = 1;
        for (i = 0; i < k; ++i) {
            
            newmui = mult((double)1/cluster_sizes[i], cluster_sums+i);
            finished &= sqr_dist(mu + i, newmui) < epsilon * epsilon;
            copydata(mu+i, newmui);
            freeVector(newmui);
        }
        finished |= iter >= maxiter;
    }
    
    
    free(cluster_sizes);
    freeVectors(cluster_sums, k);

    return mu;
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

vector* read(char* fname) {
    int d;
    int i;
    int j;
    vector* vectors;
    FILE *ifp;
    numVectors = countLines(fname);
    vectors = (vector*) calloc(numVectors, sizeof(vector));
    if (vectors == NULL) error();
    d = vectorSize(fname);
    i = 0;
    ifp = fopen(fname, "r");
    if (ifp == NULL) error();
    while(i < numVectors) {
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

void write(char* fname, vector* mu, int k) {
    FILE *fop = fopen(fname, "w");
    int i;
    int j;
    for (i = 0; i < k; ++i) {
        for (j = 0; j < mu[i].size-1; ++j) {
            fprintf(fop, "%.4f,", mu[i].data[j]);
        }
        fprintf(fop, "%.4f\n", mu[i].data[mu[i].size-1]);
    }
    fclose(fop);
}



 

int is_valid_int(char* str) {
    char* ptr = str;
    while(*ptr != '\0') {
        if (*ptr < '0' || *ptr > '9') return 0;
        ++ptr;
    }
    
    return 1;
}


int main(int argc, char* argv[]) {
    int k;
    int has_maxiter;
    int maxiter = 200;
    char* infile;
    char* outfile;
    vector* mu;
    vector* vectors;
    if (argc < 4 || argc > 5) {
        input_error();
    }
    has_maxiter = argc == 5;
    
    if (!is_valid_int(argv[1])) input_error();
    
    k = atoi(argv[1]);
    if (k <= 1) input_error();
    infile = argv[2+has_maxiter];
    outfile = argv[3+has_maxiter];
    if (has_maxiter && !is_valid_int(argv[2]) ) input_error();
    else if(has_maxiter) maxiter = atoi(argv[2]);

    if (maxiter <= 0) input_error();
    vectors = read(infile);
    if (k > numVectors) input_error();
    mu = kmeans(vectors, k, maxiter);
    write(outfile, mu, k);
    freeVectors(vectors, numVectors);
    freeVectors(mu, k);
    return 0;
}