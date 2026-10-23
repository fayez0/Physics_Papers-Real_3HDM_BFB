"""Save reproducible exact inputs; distinguish independent reference boundaries.

Six sharp reference landscapes are transformed by exact coefficient scaling
and/or doublet permutation before shifting. Their 30 targets are NOT counted
as 30 independent minimization landscapes. No reviewer-reported data is used.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json, random, hashlib, sys

ROOT = Path(__file__).resolve().parents[1]; SEED = 2026101606

def encode(p): return {k:[str(F(x)) for x in p[k]] for k in ('a','b','c','d')}
def shift(p,t): return {k:[F(x)+2*t if k in ('a','b') else F(x) for x in p[k]] for k in p}
def scale(p,t): return {k:[F(x)*t for x in p[k]] for k in p}

def build():
    rng = random.Random(SEED); rows=[]
    cohorts = [('generic',96),('vanishing_phase',48),('singular_overlap',48),('singular_radial',24),('scale_separated',24)]
    for cohort,count in cohorts:
        for j in range(count):
            p={k:[F(rng.randint(-30,30),rng.choice([5,7,9,11])) for _ in range(3)] for k in ('a','b','c','d')}
            if cohort=='vanishing_phase':
                for e in range(1+j%3): p['d'][e]=F(0)
            elif cohort=='singular_overlap':
                for e in range(1+j%3): p['c'][e]=p['d'][e]*((-1)**(j+e))
            elif cohort=='singular_radial': p['a']=[F(2),F(5),F(10)]; p['b']=[F(3),F(5),F(4)]
            elif cohort=='scale_separated':
                for k in p:
                    for e in range(3):
                        power=rng.randint(-6,6);p[k][e]*=F(10)**power
            rows.append(dict(id=f'{cohort}_{j:03d}',cohort=cohort,landscape_id=f'{cohort}_{j:03d}',couplings=encode(p),reference_kind='no_independent_global_minimum',numeric_seed=rng.randrange(2**32)))
    # These minima have separate sharp proofs in Sec. VII of the manuscript.
    singular=dict(a=[2,5,10],b=[3,5,4],c=[-2,1,4],d=[F(1,3),F(2,5),F(3,7)])
    charged=dict(a=[3]*3,b=[-2]*3,c=[1]*3,d=[F(1,2)]*3)
    rankflat=dict(a=[1]*3,b=[-1]*3,c=[2]*3,d=[0]*3)
    base=[('singularA',singular,F(43,51)),('charged_complex',charged,-F(7,60)),('rankflat',rankflat,F(0))]
    for name,p,m in base:
        for factor in [F(1),F(7,3)]:
            pp=scale(p,factor); mm=m*factor; landscape=f'{name}_scale_{str(factor).replace("/","_")}'
            for t in [-F(1,50),F(0),F(1,50),-F(1,10**20),F(1,10**20)]:
                rows.append(dict(id=f'reference_{len(rows):03d}',cohort='sharp_reference_shift',landscape_id=landscape,couplings=encode(shift(pp,t-mm)),reference_kind='independent_sharp_analytic_minimum',reference_minimum=str(t),reference_family=name,scale=str(factor),numeric_seed=rng.randrange(2**32)))
    # Deliberately charged-favoring cases; centering with the proposed solver is
    # explicitly a covariance stress test, not an independent sharp boundary.
    sys.path.insert(0,str(ROOT/'code'))
    from full_radial_infinitesimal import run_radial
    for j in range(32):
        p=dict(a=[F(rng.randint(9,20)) for _ in range(3)],b=[F(rng.randint(-3,3),10) for _ in range(3)],c=[F(rng.randint(8,16),10) for _ in range(3)],d=[F(rng.choice([-1,1])*rng.randint(5,25),100) for _ in range(3)])
        result=run_radial(**p);m=result['minimum'];landscape=f'charged_seed_{j:03d}'
        rows.append(dict(id=landscape,cohort='charged_favoring_base',landscape_id=landscape,couplings=encode(p),reference_kind='no_independent_global_minimum',numeric_seed=rng.randrange(2**32)))
        for k,t in enumerate([-F(1,50),F(0),F(1,50)]):
            rows.append(dict(id=f'charged_shift_{j:03d}_{k}',cohort='solver_centered_shift',landscape_id=landscape,couplings=encode(shift(p,t-m)),reference_kind='solver_centered_covariance_only',predicted_shifted_minimum=str(t),donor_minimum=str(m),numeric_seed=rng.randrange(2**32)))
    return rows

if __name__=='__main__':
    folder=ROOT/'validation/validation_scripts'; folder.mkdir(parents=True,exist_ok=True)
    text=json.dumps(build(),indent=2)+'\n'; (folder/'inputs.json').write_text(text)
    (folder/'input_provenance.json').write_text(json.dumps(dict(seed=SEED,rows=398,random_rows=240,charged_base_rows=32,solver_centered_shift_rows=96,reference_shift_rows=30,reference_families=3,scaled_reference_variants=6,shifted_rows_independent_landscapes=False,reviewer_data_used=False,sha256=hashlib.sha256(text.encode()).hexdigest()),indent=2)+'\n')
    print('Wrote 398 exact inputs; 30 analytic-reference shifts and 96 explicitly solver-centered covariance tests.')
