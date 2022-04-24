     
#include <math.h> 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vector.h"


PyObject* fit(PyObject* self, PyObject* args);
PyObject* PylistFromVectors(vector* vectors, int k);
vector* vectorsFromPyList(PyObject* list);

