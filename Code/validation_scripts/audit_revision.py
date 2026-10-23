"""Targeted checks of the added dispatcher, independent gradient and phase map."""
from __future__ import annotations
import itertools,json,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from exact_dispatch import classify
from field_objective import value_gradient
ROOT=Path(__file__).resolve().parents[1]

def require(value,message):
    if not value:raise ArithmeticError(message)

def run():
    rng=np.random.default_rng(2026101699);errs=[]
    for _ in range(24):
        p=rng.normal(size=(4,3));y=rng.normal(size=12);y/=np.linalg.norm(y)
        _,g=value_gradient(y,p);numeric=approx_derivative(lambda z: value_gradient(z,p)[0],y,method='3-point')
        err=float(np.max(np.abs(g-numeric)));errs.append(err);require(err<2e-8,'Direct gradient finite-difference check failed.')
    rephasings=0
    for signs in itertools.product((-1,0,1),repeat=3):
        d=np.asarray(signs,float)*np.array([.13,.19,.23]);p=np.vstack(([3,4,5],[-.4,.2,-.1],[.9,1.2,1.5],d));prod=int(np.prod(signs));target=np.abs(d);target[1]*=-1 if prod<0 else 1
        choice=None
        for t2,t3 in itertools.product((-1,1),repeat=2):
            t=np.array([1,t2,t3]);transformed=d*np.array([t[0]*t[1],t[1]*t[2],t[2]*t[0]])
            if np.array_equal(transformed,target):choice=t;break
        require(choice is not None,'No allowed discrete phase mapping.')
        z=rng.normal(size=(2,3))+1j*rng.normal(size=(2,3));zp=z*np.where(choice==1,1,1j)
        y=np.r_[z.real.ravel(),z.imag.ravel()];yp=np.r_[zp.real.ravel(),zp.imag.ravel()];pp=p.copy();pp[3]=target
        require(abs(value_gradient(y,p)[0]-value_gradient(yp,pp)[0])<2e-14,'External-convention phase map failed.');rephasings+=1
    examples=json.loads((ROOT/'examples/expected_minima.json').read_text());done=[]
    # Original example directory is a named set; values are taken directly from file.
    for name in ('charged_complex','singular_radial_positive','charged_flat','charged_nearflat_positive','charged_nearflat_negative'):
        inp=json.loads((ROOT/'examples'/f'{name}.json').read_text());result=classify(inp)
        expected={'charged_complex':-F(7,60),'singular_radial_positive':F(43,51),'charged_flat':F(0),'charged_nearflat_positive':F(1,10**20),'charged_nearflat_negative':-F(1,10**20)}[name]
        require(result['minimum']==expected,'Analytic example mismatch.');done.append(dict(id=name,minimum=str(result['minimum']),route=result['evaluation_route']))
    template=json.loads((ROOT/'examples/charged_complex.json').read_text());invalid=[]
    for name,mutate in [('binary_float',lambda p:p['a'].__setitem__(0,3.0)),('bad_rational',lambda p:p['d'].__setitem__(0,'one-half')),('missing_key',lambda p:p.pop('d')),('extra_key',lambda p:p.update(unexpected=1)),('wrong_length',lambda p:p['a'].append(2))]:
        p=json.loads(json.dumps(template));mutate(p)
        try:classify(p)
        except (TypeError,ValueError,KeyError):invalid.append(name)
        else:raise ArithmeticError(f'Invalid input accepted: {name}')
    return dict(status='passed',gradient_cases=len(errs),maximum_finite_difference_residual=max(errs),phase_sign_and_zero_patterns=rephasings,examples=done,rejected_inputs=invalid,python_optimization=not __debug__)

if __name__=='__main__':
    result=run();out=ROOT/'validation/validation_scripts'/('added_code_checks_optimized.json'if not __debug__ else'added_code_checks.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
