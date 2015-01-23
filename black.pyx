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
    """ S,K,t,v,rf,cp"""
    return tv(s,k,t,v,rf,cp)

def delta_py(double s, double k, double t,double v, double rf, double cp):
    """ S,K,t,v,rf,cp"""
    return delta(s,k,t,v,rf,cp)

def vega_py(double s, double k, double t,double v, double rf, double cp):
    """ S,K,t,v,rf,cp"""
    return vega(s,k,t,v,rf,cp)

cdef double approximate_abs_delta(double s, double k, double t,double v, double rf, double cp):
    if cp==0: ##future
        return 1
    return abs(delta(s,k,t,v,rf,cp))

@cython.cdivision(True)
@cython.boundscheck(False)
def implied_vol(double underlying, double price, double strike, double t, double rf, double cp):
    """ underlying, price, strike, tte, rf, cp """
    cdef long i = 0
    cdef double prices_guess, vol_guess = 1
    cdef double diff, delt
    
    for i in range(0,20):
        price_guess = tv(underlying,strike,t,vol_guess,rf,cp)
        diff = price - price_guess
        if abs(diff) < .001:
            return vol_guess
        vegalol = vega(underlying,strike,t,vol_guess,rf,cp)
        if vegalol<.01:
            return -1
        vol_guess += diff / vegalol
    return -1

@cython.cdivision(True)
@cython.boundscheck(False)
def implied_fut(double guess, double price, double strike, double t, double rf, double sigma, double cp):
    """ guess, price, strike, tte, rf, vol, cp """
    cdef long i
    cdef double prices_guess, underlying_guess = guess
    cdef double diff, delt
    
    if price <= .01:
        return np.NaN
    
    for i in range(20):
        price_guess = tv(underlying_guess,strike,t,sigma,rf,cp)
        diff = price - price_guess
        if abs(diff) < .0001:
            return long(underlying_guess*10000) / 10000.0
        delt = delta(underlying_guess,strike,t,sigma,rf,cp)
        underlying_guess += diff / delt
    return np.NaN

#quick functions to go back and from call tvs to vols and vv.
def vols_to_tvs(spot,ks,tte,vs,ir=.02,cp=1):
    """ spot, ks, tte, vs, ir, cp """
    if len(ks)!=len(vs):
        raise "StrikeVols Unaligned error"
    return [tv(spot,p[0],tte,p[1],ir,cp)for p in zip(ks,vs)]

def tvs_to_vols(spot,ks,tte,tvs,ir=.02,cp=1):
    """ spot, ks, tte, tvs, ir, cp """
    if len(ks)!=len(tvs):
        raise "StrikeVols Unaligned error"
    return [implied_vol(spot,p[0],p[1],tte,ir,cp)for p in zip(tvs,ks)]

