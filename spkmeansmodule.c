#define PY_SSIZE_T_CLEAN  
#include <Python.h> 
#include "kmeans.h"
#include "spkmeans.h"


PyObject* PyListFromMatrix(double** matrix, int n) {
    
    PyObject* list;
    int i, j;

    list = PyList_New(n);
    for (i = 0; i < n; ++i) {
        PyList_SetItem(list, i, PyList_New(n));
        for (j = 0; j < n; ++j) {
            PyList_SetItem(PyList_GetItem(list, i), j, PyFloat_FromDouble(matrix[i][j]));
        }
    }
    return list;
}

double** matrixFromPyList(PyObject* list) {
    Py_ssize_t size;
    double** matrix;
    PyObject* row;
    size = PyList_Size(list);
    
    int i, j;
    matrix = new_matrix(size);

    for (i = 0; i < size; ++i) {
        row = PyList_GetItem(list, i);
    
        for (j = 0; j < size; ++j) {
            matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }
    return matrix;
}

PyObject* PyListFromEigen(eigen eigens, int n) {
    
    PyObject* ret = PyList_New(2);
    int i;
    PyList_SetItem(ret, 0, PyListFromMatrix(eigens.vecs, n));
    PyList_SetItem(ret, 1 , PyList_New(n));
    for (i = 0; i < n; ++i) {
        PyList_SetItem(PyList_GetItem(ret, 1), i, PyFloat_FromDouble(eigens.vals[i]));
    }
    return ret;
}



PyObject* wam(PyObject* self, PyObject* args) {
    PyObject* py_data;
    PyObject* returnList;
    vector* vectors;
    double** matrix;
    if(!PyArg_ParseTuple(args, "O", &py_data)) return NULL;
    vectors = vectorsFromPyList(py_data);

    matrix = wam_c(vectors, PyList_Size(py_data));
    
    returnList = PyListFromMatrix(matrix, PyList_Size(py_data));
    free_matrix(matrix);
    freeVectors(vectors, PyList_Size(py_data));
    return returnList;
}

PyObject* ddg(PyObject* self, PyObject* args) {
    PyObject* py_data;
    PyObject* returnList;
    vector* vectors;
    double** matrix;
    if(!PyArg_ParseTuple(args, "O", &py_data)) return NULL;
    vectors = vectorsFromPyList(py_data);

    matrix = ddg_c(vectors, PyList_Size(py_data));
    
    returnList = PyListFromMatrix(matrix, PyList_Size(py_data));
    free_matrix(matrix);
    freeVectors(vectors, PyList_Size(py_data));
    return returnList;
}

PyObject* lnorm(PyObject* self, PyObject* args) {
    PyObject* py_data;
    PyObject* returnList;
    vector* vectors;
    double** matrix;
    if(!PyArg_ParseTuple(args, "O", &py_data)) return NULL;
    vectors = vectorsFromPyList(py_data);

    matrix = lnorm_c(vectors, PyList_Size(py_data));
    
    returnList = PyListFromMatrix(matrix, PyList_Size(py_data));
    free_matrix(matrix);
    freeVectors(vectors, PyList_Size(py_data));
    return returnList;
}

PyObject* jacobi(PyObject* self, PyObject* args) {
    PyObject* py_data;
    PyObject* returnList;
    double** matrix;
    eigen ret;
    if(!PyArg_ParseTuple(args, "O", &py_data)) return NULL;
    matrix = matrixFromPyList(py_data);

    ret = jacobi_c(matrix, PyList_Size(py_data));

    returnList = PyListFromEigen(ret, PyList_Size(py_data));
    free_matrix(ret.vecs);
    free(ret.vals);
    
    return returnList;
}

static PyMethodDef spkmeansMethods[] = {
    {"ddg",
        (PyCFunction) ddg,
        METH_VARARGS,
        PyDoc_STR("Calculate the diagonal matrix")},
    {"wam",
        (PyCFunction) wam,
        METH_VARARGS,
        PyDoc_STR("Calculate the weighted matrix")},
    {"lnorm",
        (PyCFunction) lnorm,
        METH_VARARGS,
        PyDoc_STR("Calculate the lnorm matrix")},
    {"jacobi",
        (PyCFunction) jacobi,
        METH_VARARGS,
        PyDoc_STR("Find eigen vectors and eigenvalues")},
    {"fit",
        (PyCFunction) fit,
        METH_VARARGS,
        PyDoc_STR("Find eigen vectors and eigenvalues")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL,
    -1,
    spkmeansMethods
};

PyMODINIT_FUNC 
PyInit_spkmeans(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {return NULL;}
    return m;
}

