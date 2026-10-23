"""Reproducible validation; 40-start local searches are NOT global certificates.

Timing is in-process after imports. Exact always-deformed and regular-first
routes run in the same process. Numerical optimization runs only against the
independent original-field module, using coefficient normalization for scale.
Every attempt and failure flag is saved; no failed start is discarded silently.
"""
from __future__ import annotations
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[name]='1'
import argparse, concurrent.futures, hashlib, importlib.metadata, json, platform, sys, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from exact_dispatch import classify, run_radial
from field_objective import value_gradient
ROOT=Path(__file__).resolve().parents[1]


def work(task):
    row,starts,maxiter=task; p=row['couplings']; t=time.perf_counter()
    try:
        exact=run_radial(**p,deform=True); full_time=time.perf_counter()-t; fast=classify(p)
        if exact['minimum'] != fast['minimum']: raise ArithmeticError('Exact dispatch disagreement.')
        if 'reference_minimum' in row and exact['minimum'] != F(row['reference_minimum']): raise ArithmeticError('Independent analytic reference disagreement.')
        coefficients=np.array([[float(F(x)) for x in p[k]] for k in ('a','b','c','d')]); scale=max(1.,float(np.max(np.abs(coefficients)))); pn=coefficients/scale
        rng=np.random.default_rng(row['numeric_seed']); numerical=[]; t=time.perf_counter(); best_y=None; best=float('inf')
        for j in range(starts):
            y=rng.normal(size=12); y/=np.linalg.norm(y)
            try:
                sol=minimize(value_gradient,y,args=(pn,),jac=True,method='BFGS',options={'gtol':1e-10,'maxiter':maxiter})
                val,grad=value_gradient(sol.x,pn); yn=sol.x/np.linalg.norm(sol.x)
                rec=dict(start=j,value=val*scale,scaled_value=val,success=bool(sol.success),status=int(sol.status),message=str(sol.message),iterations=int(sol.nit),gradient_norm=float(np.linalg.norm(grad)),norm=float(np.dot(yn,yn)))
                if np.isfinite(val) and val<best:best=val;best_y=yn.tolist()
            except Exception as exc:rec=dict(start=j,success=False,status='exception',message=f'{type(exc).__name__}: {exc}')
            numerical.append(rec)
        numeric_time=time.perf_counter()-t
        if not np.isfinite(best):raise ArithmeticError('All numerical starts failed to produce a finite value.')
        residual=best-float(exact['minimum'])/scale
        return dict(**row,completed=True,minimum=str(exact['minimum']),classification=exact['classification'],branch=exact['best']['label'],exact_deformed_seconds=full_time,
                    exact_dispatch_seconds=fast['elapsed_seconds'],dispatch_route=fast['evaluation_route'],fallback_reason=fast['fallback_reason'],
                    numeric_seconds=numeric_time,coefficient_scale=scale,best_numeric=best*scale,best_fields=best_y,
                    difference_scaled=residual,difference_absolute=residual*scale,scaled_comparison_tolerance=1e-7,
                    comparison='agreement_at_stated_tolerance' if abs(residual)<=1e-7 else ('lower_point_found' if residual<0 else 'local_search_above_exact_result'),starts=numerical)
    except Exception as exc:return dict(**row,completed=False,error=f'{type(exc).__name__}: {exc}')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',type=Path,default=ROOT/'validation/validation_scripts/inputs.json');parser.add_argument('--output',type=Path,default=ROOT/'validation/validation_scripts/records.jsonl')
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--starts',type=int,default=40);parser.add_argument('--maxiter',type=int,default=1500);parser.add_argument('--start',type=int,default=0);parser.add_argument('--stop',type=int);args=parser.parse_args()
    rows=json.loads(args.inputs.read_text())[args.start:args.stop];args.output.parent.mkdir(parents=True,exist_ok=True)
    env=dict(python=sys.version,platform=platform.platform(),processor=platform.processor(),cpu_count=os.cpu_count(),workers=args.workers,starts_per_input=args.starts,maxiter=args.maxiter,
             method='BFGS',gtol=1e-10,scaled_comparison_tolerance=1e-7,timing='in-process; imports excluded; concurrent worker load',seed_source='per-input numeric_seed',
             proof_status='Numerical agreement is not a global certificate; only separately identified sharp references are independent exact minima.',input_sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest())
    for name in ('numpy','scipy','sympy','mpmath'):env[name]=importlib.metadata.version(name)
    (args.output.parent/'environment.json').write_text(json.dumps(env,indent=2)+'\n');t=time.perf_counter();fail=0
    with args.output.open('w') as stream, concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,r in enumerate(pool.map(work,[(r,args.starts,args.maxiter) for r in rows])):
            stream.write(json.dumps(r,default=str)+'\n');stream.flush();fail+=int(not r['completed'] or r.get('comparison')=='lower_point_found')
            if i%10==0:print(json.dumps(dict(done=i+1,id=r['id'],completed=r['completed'],comparison=r.get('comparison'),elapsed_seconds=time.perf_counter()-t)),flush=True)
    print(json.dumps(dict(rows=len(rows),failures=fail,total_wall_seconds=time.perf_counter()-t)),flush=True)
    return int(fail>0)
if __name__=='__main__':raise SystemExit(main())
