#include "python_wrapper.h"
#include "fpotencia.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/************************************************************************************************
Pure C++ function implementation
 ************************************************************************************************/

/*
Returns the probability of x given sigma[s] and the mean [xm]
 */
static double normal1(double x, double xm, double s) {

    Load ld23("Load23", 23, 8.7, 6.7);


    double B[] = {1.330274429, -1.821255978, 1.781477937, -0.356563782, 0.319381530};

    double pp = 0.2316419;
    double y = (x - xm) / s; // reduced variable
    double zy = exp(-y * y / 2.0) / 2.506628275;

    //calculate qx by Horner's method
    double t = 1.0 / (1.0 + pp * y);
    double po = 0.0;
    for (double Bi : B)
        po = po * t + Bi;
    po = po * t;
    double qx = zy * po; //1 - probability = tail
    //double px = 1.0 - qx; //probability
    return 1.0 - qx; //
}

/*
Returns the two tails of the normal law centered in 0
 */
static double normal2(double x, double s) {
    double xm = 0.;
    double B[] = {1.330274429, -1.821255978, 1.781477937, -0.356563782, 0.319381530};

    double pp = 0.2316419;
    double y = (x - xm) / s; // reduced variable
    double zy = exp(-y * y / 2.0) / 2.506628275;

    //calculate qx by Horner's method
    double t = 1.0 / (1.0 + pp * y);
    double po = 0.0;
    for (double Bi : B)
        po = po * t + Bi;
    po = po * t;
    double qx = zy * po; //1 - probability = tail
    //*px = 1.0 - *qx; #probability
    return 2.0 * qx;
}

/*
Distance between 2, 2D points.
px: x coordinate of the point 'p'
py: y coordinate of the point 'p'
sx: x coordinate of the point 's'
sy: y coordinate of the point 's'
 */
double euclidean_distance(double px, double sx, double py, double sy) {
    return sqrt((px - sx)*(px - sx) + (py - sy)*(py - sy));
}

/************************************************************************************************
Python wrapper functions (where possible)
 ************************************************************************************************/

static PyObject *load(PyObject *self, PyObject *args) {


    Load ld23("Load23", 23, 8.7, 6.7);

    return PyObject_FROM_GC(ld23);
}

/*
Example: trace of a matrix
 */
static PyObject *trace(PyObject *self, PyObject *args) {
    PyObject *input;
    PyArrayObject *array;
    double sum;
    int i, n;

    if (!PyArg_ParseTuple(args, "O", &input))
        return NULL;

    array = (PyArrayObject *) PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 2, 2);

    if (array == NULL)
        return NULL;

    n = array->dimensions[0];

    if (n > array->dimensions[1])
        n = array->dimensions[1];

    sum = 0.;

    for (i = 0; i < n; i++)
        sum += *(double *) (array->data + i * array->strides[0] + i * array->strides[1]);

    Py_DECREF(array);

    return PyFloat_FromDouble(sum);
}


/************************************************************************************************
Python module construction
 ************************************************************************************************/

//Array defining the visible functions
PyMethodDef pymethods[] = {
    { "trace", (PyCFunction) trace, METH_VARARGS, NULL},
    { NULL, NULL, 0, NULL} /*Always present*/
};


//constructor

PyMODINIT_FUNC initfPotencia(void) {
    PyObject *m;

    m = Py_InitModule("fPotencia", pymethods);

    import_array(); //use when using 
    if (m == NULL)
        return;

}