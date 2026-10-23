"""Export all 28 stationary branches, including rejections, over Q(epsilon).

The certificate requires the full-domain branch-completeness theorem. It is
not an unconditional sum-of-squares certificate. The verifier is a separate
module, imports no solver, and never solves or minimizes anything.
"""
from fractions import Fraction as F
import full_infinitesimal as f

def pack(x):
    def coeff(poly):
        terms=poly.to_dict(); degree=max((k[0] for k in terms),default=0)
        return [str(terms.get((i,),0)) for i in range(degree+1)]
    return {'numerator':coeff(x.numer),'denominator':coeff(x.denom)}

def export_decision(a,b,c,d):
    raw=dict(a=a,b=b,c=c,d=d)
    if any(len(v)!=3 for v in raw.values()): raise ValueError('Expected four length-three arrays')
    if any(isinstance(x,(float,bool)) for v in raw.values() for x in v): raise TypeError('Exact rational inputs required')
    raw={k:list(map(F,v)) for k,v in raw.items()}
    p={k:[f.scalar(x)+f.EPS*f.PERTURBATION[k][i] for i,x in enumerate(v)] for k,v in raw.items()}
    a,b,c,d=[p[k] for k in ('a','b','c','d')]; rows=[]; limits=[]
    A=[[f.Z]*3 for _ in range(3)]
    for i in range(3):A[i][i]=a[i]
    for e,(i,j) in enumerate(f.EDGES):A[i][j]=A[j][i]=b[e]
    R=f.inverse(A); w=f.mv(R,[f.O]*3)
    alpha=[c[e]+d[e] for e in range(3)]; beta=[c[e]-d[e] for e in range(3)]
    for sector in ('neutral','charged'):
        if sector=='neutral':diag=[x/2 for x in a];off=[(b[e]+c[e])/2 for e in range(3)];delta=[x/2 for x in d];weight=[-f.O]*3
        else:diag=[R[i][i] for i in range(3)];off=[R[i][j]+(1/alpha[e]+1/beta[e])/2 for e,(i,j) in enumerate(f.EDGES)];delta=[(1/alpha[e]-1/beta[e])/2 for e in range(3)];weight=w
        for support,Q,kind,signs in f.branches(diag,off,delta):
            z=f.solve(f.bordered(Q,support,weight),[f.Z]*len(support)+[f.O]);s=[f.Z]*3
            for i,v in zip(support,z):s[i]=v
            scalar=z[-1];row=dict(sector=sector,support=list(support),kind=kind,signs=signs,s=[pack(v) for v in s],stationary_scalar=pack(scalar));reason=None;value=None
            if any(f.sgn(s[i])<=0 for i in support):reason='nonpositive_support'
            if reason is None:
                cos=f.phase_cosines(s,kind,signs,delta)
                if cos is None:reason='unrealizable_triangle'
            if reason is None and sector=='neutral':value=scalar
            if reason is None and sector=='charged':
                T=f.dot(w,s)+scalar*sum(w,f.Z)
                if not T:reason='zero_normalization'
                else:
                    re=f.mv(R,s);trE2=f.dot(re,re)+sum((s[i]*s[j]*((1/alpha[e]**2+1/beta[e]**2)+(1/alpha[e]**2-1/beta[e]**2)*cos[e]) for e,(i,j) in enumerate(f.EDGES)),f.Z)
                    trK2=(trE2+2*scalar*f.dot(re,w)+scalar**2*f.dot(w,w))/T**2
                    if f.sgn(1-trK2)<0:reason='indefinite_gram'
                    else:value=scalar/(2*T)
            row.update(accepted=reason is None,rejection=reason)
            if value is not None:row.update(value=pack(value),value_limit=str(f.limit_zero(value)));limits.append(f.limit_zero(value))
            rows.append(row)
    if len(rows)!=28 or not limits:raise RuntimeError('Incomplete certificate construction')
    m=min(limits)
    return dict(schema='weinberg-full-decision-v1',theorem_dependency='FULL_DOMAIN_PROOF.md Sections 2-7',coefficients={k:list(map(str,v)) for k,v in raw.items()},regularizer=f.PERTURBATION,rows=rows,minimum=str(m),classification='positive' if m>0 else 'flat' if m==0 else 'negative')
