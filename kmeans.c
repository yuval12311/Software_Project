#define PY_SSIZE_T_CLEAN  
#include <Python.h>      
#include <math.h> 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
vector* kmeans(vector* vectors, vector* mu, int k, double epsilon, int maxiter) {
    /* vector* mu = (vector*) calloc(k, sizeof(vector)); */
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
    /* for(bob = 0; bob < k; ++bob) {
        (mu+bob)->data =(double*) malloc(sizeof(double)*vectors->size);
        if ((mu+bob)->data == NULL) error(); 
        copydata(mu+bob, vectors+bob);
    } */
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

vector* vectorsFromPyList(PyObject* list) {
    Py_ssize_t size;
    vector* vectors;
    PyObject* row;
    size = PyList_Size(list);
    
    int i, j;
    vectors =(vector*) malloc(size * sizeof(vector));
    for (i = 0; i < size; ++i) {
        row = PyList_GetItem(list, i);
        vectors[i].size = PyList_Size(row);
        vectors[i].data = malloc(vectors[i].size * sizeof(double));
        for (j = 0; j < vectors[i].size; ++j) {
            vectors[i].data[j] = PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }
    return vectors;
    /* maybe nned to decref */
}


PyObject* PylistFromVectors(vector* vectors, int k) {
    PyObject* list;
    int i, j;
    int d = vectors->size;

    list = PyList_New(k);
    for (i = 0; i < k; ++i) {
        PyList_SetItem(list, i, PyList_New(d));
        for (j = 0; j < d; ++j) {
            PyList_SetItem(PyList_GetItem(list, i), j, PyFloat_FromDouble(vectors[i].data[j]));
        }
    }
    return list;
}

PyObject* fit(PyObject* self, PyObject* args) {
    int k;
    double epsilon;
    int maxiter;
    PyObject* py_initial_mu;
    PyObject* py_data;
    PyObject* returnList;
    vector* mu;
    vector* data;
    if(!PyArg_ParseTuple(args, "idiOO", &k, &epsilon, &maxiter, &py_initial_mu, &py_data)) return NULL;
    mu = vectorsFromPyList(py_initial_mu);
    data = vectorsFromPyList(py_data);
    numVectors = PyList_Size(py_data);

    kmeans(data, mu, k, epsilon, maxiter);

    returnList = PylistFromVectors(mu, k);

    freeVectors(mu, k);
    freeVectors(data, numVectors);




    return returnList;

}

static PyMethodDef kmeansMethods[] = {
    {"fit",
        (PyCFunction) fit,
        METH_VARARGS,
        PyDoc_STR("Runs the Kmeans algorithm on a dataset with given k and initial centroids")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "kmeans",
    NULL,
    -1,
    kmeansMethods
};

PyMODINIT_FUNC 
PyInit_kmeans(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {return NULL;}
    return m;
}






 
