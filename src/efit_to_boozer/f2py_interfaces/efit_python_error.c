#include <Python.h>

/* f2py holds the GIL and checks PyErr_Occurred after the Fortran call. */
void libneo_efit_python_error(const char *message)
{
    PyErr_SetString(PyExc_ValueError, message);
}
