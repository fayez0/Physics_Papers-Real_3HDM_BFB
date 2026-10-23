"""Direct 12-real-field objective and gradient; no Gram/branch-solver imports.

Columns of u + i*v are the three doublets. Coefficients are floats here solely
for numerical cross-checks. The returned value is V4 / (sum norm^2)^2.
This independent implementation does not certify global optimality.
"""
from __future__ import annotations
import numpy as np

PAIRS = ((0,1),(1,2),(2,0))

def value_gradient(y: np.ndarray, coefficients: np.ndarray) -> tuple[float, np.ndarray]:
    y = np.asarray(y, dtype=float); p = np.asarray(coefficients, dtype=float)
    if y.shape != (12,) or p.shape != (4,3): raise ValueError('Need 12 fields and four coefficient triples.')
    u, v = y[:6].reshape(2,3), y[6:].reshape(2,3)
    r = np.sum(u*u+v*v,axis=0); S = float(r.sum())
    if not np.isfinite(S) or S <= 0: raise ValueError('Field norm must be finite and nonzero.')
    a,b,c,d = p; V = float(np.dot(a,r*r)/2); dr = a*r
    for e,(i,j) in enumerate(PAIRS): V += b[e]*r[i]*r[j]; dr[i] += b[e]*r[j]; dr[j] += b[e]*r[i]
    gu,gv = 2*u*dr, 2*v*dr
    for e,(i,j) in enumerate(PAIRS):
        real = float(np.dot(u[:,i],u[:,j])+np.dot(v[:,i],v[:,j]))
        imag = float(np.dot(u[:,i],v[:,j])-np.dot(v[:,i],u[:,j]))
        cr,ci = c[e]+d[e],c[e]-d[e]; V += cr*real*real+ci*imag*imag
        gu[:,i] += 2*cr*real*u[:,j]+2*ci*imag*v[:,j]; gv[:,i] += 2*cr*real*v[:,j]-2*ci*imag*u[:,j]
        gu[:,j] += 2*cr*real*u[:,i]-2*ci*imag*v[:,i]; gv[:,j] += 2*cr*real*v[:,i]+2*ci*imag*u[:,i]
    grad = np.r_[gu.ravel(),gv.ravel()]/S**2-4*V*y/S**3
    return V/S**2, grad
