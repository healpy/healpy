#ifndef HEALPY_UTILS_H
#define HEALPY_UTILS_H

#define healpyAssertType(testval,msg) \
do { if (testval); else { PyErr_SetString(PyExc_TypeError, msg); return NULL; } } while(0)
#define healpyAssertValue(testval,msg) \
do { if (testval); else { PyErr_SetString(PyExc_ValueError, msg); return NULL; } } while(0)

#define CP_(str)(char *)(str)

#endif
