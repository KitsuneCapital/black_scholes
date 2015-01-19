import numpy as np
cimport numpy as np
import pandas as pd
from libc.math cimport exp, sqrt, pow, log, erf, abs, M_PI
cimport cython

ctypedef np.double_t DTYPE_t


cdef extern from "black.h":
    double tv(double s, double k, double t,double v, double rf, double cp)
    double vega(double s, double k, double t,double v, double rf, double cp)
    double delta(double s, double k, double t,double v, double rf, double cp)

def tv_py(double s, double k, double t,double v, double rf, double cp):
    return tv(s,k,t,v,rf,cp)

def delta_py(double s, double k, double t,double v, double rf, double cp):
    return delta(s,k,t,v,rf,cp)

def vega_py(double s, double k, double t,double v, double rf, double cp):
    return vega(s,k,t,v,rf,cp)

