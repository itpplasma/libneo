/* fpm compiles every source; without Python headers this unit is empty. */
#if defined(__has_include)
#if __has_include(<Python.h>)
#define LIBNEO_HAVE_PYTHON_H 1
#endif
#endif

#ifdef LIBNEO_HAVE_PYTHON_H
#include <Python.h>

/* f2py holds the GIL and checks PyErr_Occurred after the Fortran call. */
void libneo_efit_python_error(const char *message)
{
    PyErr_SetString(PyExc_ValueError, message);
}
#endif
