"""Independent field-space numerical cross-check; not a proof or certificate."""
import json,time
import numpy as np
from scipy.optimize import minimize
from prototype import Couplings,minimum,value,L,EDGES

def objective(y,p):
    X=y[:6].reshape(2,3)+1j*y[6:].reshape(2,3); S=float(np.vdot(X,X).real); K=X.conj().T@X; V=value(K,p)
    G=2*(X@L(K,p))/S**2-4*V*X/S**3
    return V/S**2,np.r_[G.real.ravel(),G.imag.ravel()]

def run(count=100,starts=8,seed=20261016):
    rng=np.random.default_rng(seed); records=[]; maxerr=0.; t0=time.monotonic(); charged=0; tri=0; missed=0; failures=0
    for k in range(count):
        p=Couplings(rng.integers(2,25,3)/5,rng.integers(-30,30,3)/7,rng.integers(-20,25,3)/9,rng.integers(-22,22,3)/11)
        # Integer generation can hit exact degeneracies; keep them visible and exclude from regular-case comparison.
        best,pts,unresolved=minimum(p)
        if unresolved: records.append({'index':k,'excluded':str(unresolved)}); continue
        direct=[]
        for j in range(starts):
            x=rng.normal(size=12); sol=minimize(objective,x,args=(p,),jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':700}); direct.append(float(sol.fun)); failures+=int(not sol.success)
        numeric=min(direct); err=abs(numeric-best['value']); maxerr=max(maxerr,err); charged+=int(best['sector'].startswith('charged')); tri+=int(best['sector']=='charged' and best['branch']=='triangle'); missed+=int(numeric-best['value']>1e-7)
        if numeric<best['value']-1e-7: raise AssertionError(('missing lower candidate',k,numeric,best['value'],p))
        records.append({'index':k,'a':p.a.tolist(),'b':p.b.tolist(),'c':p.c.tolist(),'d':p.d.tolist(),'exact_formula_float':best['value'],'direct_multistart':numeric,'difference':numeric-best['value'],'sector':best['sector'],'branch':best['branch'],'support':best['support'],'stationary_residual':best['stationary_residual']})
    report={'seed':seed,'points':count,'starts_per_point':starts,'max_absolute_difference':maxerr,'charged_minima':charged,'charged_triangle_minima':tri,'multistart_above_formula':missed,'optimizer_non_success_flags':failures,'elapsed_seconds':time.monotonic()-t0,'records':records}
    return report

if __name__=='__main__':
    import sys
    count=int(sys.argv[1]) if len(sys.argv)>1 else 100
    report=run(count=count); open(f'numeric_audit_{count}.json','w').write(json.dumps(report,indent=2)); print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
