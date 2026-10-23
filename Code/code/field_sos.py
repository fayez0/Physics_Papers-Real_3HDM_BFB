"""Independent original-field rational SOS lower bounds; no BFB solver imports.

The optional numerical search uses alternating projections solely to propose a
Gram matrix. Only exact coefficient matching and exact rational LDL factors
establish a certificate. Failure to find a certificate makes NO classification.
Coordinates: Phi_i=(x[4i]+i*x[4i+1], x[4i+2]+i*x[4i+3]).
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from collections import defaultdict
import numpy as np

EDGES=((0,1),(1,2),(2,0))

def frac(x):
    if isinstance(x,(float,np.floating)): raise TypeError('Exact rational inputs required')
    return F(x)

def addto(a,b,scale=F(1)):
    for mon,c in b.items():
        a[mon]=a.get(mon,F(0))+scale*c
        if not a[mon]: del a[mon]
    return a

def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            k=tuple(sorted(x+y)); out[k]=out.get(k,F(0))+c*d
    return {k:v for k,v in out.items() if v}

def field_polynomial(a,b,c,d,lower=0):
    p=[list(map(frac,v)) for v in (a,b,c,d)]
    if any(len(v)!=3 for v in p): raise ValueError('Expected four coefficient triples')
    a,b,c,d=p; r=[{(4*i+j,4*i+j):F(1) for j in range(4)} for i in range(3)]; S={}; out={}
    for i in range(3): addto(S,r[i]); addto(out,mul(r[i],r[i]),a[i]/2)
    for e,(i,j) in enumerate(EDGES):
        re={}; im={}
        for k in (0,2):
            for u,v,sg in ((0,0,1),(1,1,1)): addto(re,{tuple(sorted((4*i+k+u,4*j+k+v))):F(sg)})
            for u,v,sg in ((0,1,1),(1,0,-1)): addto(im,{tuple(sorted((4*i+k+u,4*j+k+v))):F(sg)})
        addto(out,mul(r[i],r[j]),b[e]); addto(out,mul(re,re),c[e]+d[e]); addto(out,mul(im,im),c[e]-d[e])
    addto(out,mul(S,S),-frac(lower)); return out

def bases(degree=2):
    groups=defaultdict(list)
    for ij in combinations_with_replacement(range(12),degree):
        key=tuple(sum(i//4==k for i in ij)%2 for k in range(3))+(sum(i%2 for i in ij)%2,)
        groups[key].append(ij)
    return [groups[k] for k in sorted(groups)]

def exact_ldl(matrix):
    n=len(matrix); A=[list(map(frac,row)) for row in matrix]
    if any(len(row)!=n for row in A) or any(A[i][j]!=A[j][i] for i in range(n) for j in range(n)): raise ValueError('Not a symmetric square matrix')
    L=[[F(i==j) for j in range(n)] for i in range(n)]; diag=[]
    for j in range(n):
        pivot=A[j][j]-sum((L[j][k]**2*diag[k] for k in range(j)),F(0))
        if pivot<0: raise ValueError('Negative exact LDL pivot')
        diag.append(pivot)
        for i in range(j+1,n):
            rem=A[i][j]-sum((L[i][k]*L[j][k]*diag[k] for k in range(j)),F(0))
            if not pivot:
                if rem: raise ValueError('Zero pivot with nonzero residual column')
            else: L[i][j]=rem/pivot
    return L,diag

def verify_sos(cert):
    """Exact verifier: checks every monomial and every rational PSD factor."""
    if cert.get('schema')!='3hdm-original-fields-sos-v1': raise ValueError('Unknown schema')
    p=cert['coefficients']; target=field_polynomial(**p,lower=cert['lower_bound']); got={}; ranks=[]; multiplier=cert.get('sphere_multiplier_power',0)
    if not isinstance(multiplier,int) or multiplier<0 or multiplier>4: raise ValueError('Unsupported multiplier')
    for _ in range(multiplier): target=mul(target,{(i,i):F(1) for i in range(12)})
    basis_blocks=bases(2+multiplier)
    if len(cert['blocks'])!=len(basis_blocks): raise ValueError('Missing blocks')
    for record,basis in zip(cert['blocks'],basis_blocks):
        if record['monomials']!=list(map(list,basis)): raise ValueError('Unexpected monomials')
        Q=[list(map(frac,row)) for row in record['gram']]; L,D=exact_ldl(Q); ranks.append(sum(v>0 for v in D))
        if len(Q)!=len(basis): raise ValueError('Wrong block dimensions')
        for i,u in enumerate(basis):
            for j,v in enumerate(basis): addto(got,{tuple(sorted(u+v)):Q[i][j]})
    if got!=target: raise ValueError('Exact polynomial identity failed')
    return dict(verified=True,lower_bound=str(frac(cert['lower_bound'])),original_real_variables=12,quartic_monomials=len(target),block_orders=[len(b) for b in basis_blocks],sphere_multiplier_power=multiplier,block_ranks=ranks,independent_of_branch_theorem=True)

class GramSearch:
    def __init__(self,multiplier=0):
        self.multiplier=multiplier; self.basis=bases(2+multiplier); self.sizes=[len(b) for b in self.basis]; self.ends=np.cumsum([0]+[n*n for n in self.sizes]); mons=sorted({tuple(sorted(u+v)) for B in self.basis for u in B for v in B}); self.mons=mons; self.index={m:i for i,m in enumerate(mons)}
        self.labels=np.array([self.index[tuple(sorted(u+v))] for B in self.basis for u in B for v in B]); self.count=np.bincount(self.labels,minlength=len(mons)); self.length=len(self.labels)
    def affine(self,x,b): return x+((b-np.bincount(self.labels,weights=x,minlength=len(self.mons)))/self.count)[self.labels]
    def cone(self,x,floor=0):
        z=np.empty_like(x)
        for k,n in enumerate(self.sizes):
            A=x[self.ends[k]:self.ends[k+1]].reshape(n,n); eig,U=np.linalg.eigh((A+A.T)/2); z[self.ends[k]:self.ends[k+1]]=((U*np.maximum(eig,floor))@U.T).ravel()
        return z
    def mineig(self,x): return min(np.linalg.eigvalsh(x[self.ends[k]:self.ends[k+1]].reshape(n,n))[0] for k,n in enumerate(self.sizes))
    def rationalize(self,x,poly,coefficients,lower,denominator=10**10):
        q=[F(int(round(v*denominator)),denominator) for v in x]; residual=[poly.get(m,F(0)) for m in self.mons]
        for j,k in enumerate(self.labels): residual[k]-=q[j]
        for j,k in enumerate(self.labels): q[j]+=residual[k]/int(self.count[k])
        blocks=[]
        for k,B in enumerate(self.basis):
            n=len(B); flat=q[self.ends[k]:self.ends[k+1]]; Q=[flat[i*n:(i+1)*n] for i in range(n)]; exact_ldl(Q); blocks.append(dict(monomials=list(map(list,B)),gram=[[str(v) for v in row] for row in Q]))
        cert=dict(schema='3hdm-original-fields-sos-v1',coefficients={k:list(map(str,v)) for k,v in coefficients.items()},lower_bound=str(frac(lower)),sphere_multiplier_power=self.multiplier,blocks=blocks)
        verify_sos(cert); return cert
    def search(self,p,lower,maxiter=3000,tolerance=1e-9):
        poly=field_polynomial(**p,lower=lower)
        for _ in range(self.multiplier): poly=mul(poly,{(i,i):F(1) for i in range(12)})
        if set(poly)-set(self.mons): raise ValueError('Polynomial violates declared parity symmetries')
        b=np.array([float(poly.get(m,0)) for m in self.mons]); scale=max(1.,np.max(np.abs(b))); b=b/scale; x=self.affine(np.zeros(self.length),b); log=[]; floor=1e-8
        for step in range(maxiter):
            y=self.cone(x,floor); z=self.affine(2*y-x,b); x=x+z-y
            if step%25==0 or step==maxiter-1:
                q=self.affine(y,b); eigen=self.mineig(q); log.append(dict(iteration=step,min_eigenvalue=float(eigen*scale),residual=float(np.linalg.norm(y-q)*scale)))
                if eigen>tolerance:
                    try:
                        cert=self.rationalize(q*scale,poly,p,lower); return dict(status='exact_sos_certificate',iterations=step+1,certificate=cert,log=log)
                    except ValueError: pass
        return dict(status='no_certificate_found',iterations=maxiter,log=log)

if __name__=='__main__':
    import json,sys,time
    p=dict(a=[3]*3,b=[-2]*3,c=[1]*3,d=['1/2']*3); lower=F(-7,60)-F(1,1000); t=time.monotonic(); result=GramSearch().search(p,lower); print(json.dumps({k:v for k,v in result.items() if k not in ('certificate','log')},indent=2),time.monotonic()-t); print(result['log'][-3:]);
    if 'certificate' in result: open(sys.argv[1] if len(sys.argv)>1 else 'test_sos.json','w').write(json.dumps(result['certificate'],indent=2))

def verify_separator(cert):
    """Prove the specified quartic is NOT SOS using a rational moment functional."""
    if cert.get('schema')!='3hdm-sos-separator-v1': raise ValueError('Unknown separator schema')
    B=bases(); expected={tuple(sorted(u+v)) for block in B for u in block for v in block}; rows=cert['moment_values']; moments={tuple(row['monomial']):frac(row['value']) for row in rows}
    if len(rows)!=len(moments) or set(moments)!=expected: raise ValueError('Incomplete moment functional')
    ranks=[]
    for block in B:
        matrix=[[moments[tuple(sorted(u+v))] for v in block] for u in block]; _,D=exact_ldl(matrix); ranks.append(sum(v>0 for v in D))
    p=cert['coefficients']; target=field_polynomial(**p,lower=cert['lower_bound']); pairing=sum((v*moments[m] for m,v in target.items()),F(0))
    if pairing>=0 or pairing!=frac(cert['separation']): raise ValueError('No strict exact separation')
    S2=mul({(i,i):F(1) for i in range(12)},{(i,i):F(1) for i in range(12)}); base=field_polynomial(**p); pair=lambda f:sum((v*moments[m] for m,v in f.items()),F(0)); norm=pair(S2)
    if norm<=0: raise ValueError('Invalid normalization')
    return dict(verified=True,separation=str(pairing),sos_bound_upper=str(pair(base)/norm),block_ranks=ranks,full_real_variables=12,independent_of_branch_theorem=True)

def exact_field_value(p,x):
    """Evaluate the ORIGINAL doublet expression with rational real coordinates."""
    if len(x)!=12: raise ValueError('Expected twelve real coordinates')
    x=list(map(frac,x)); a,b,c,d=[list(map(frac,p[k])) for k in ('a','b','c','d')]; r=[sum(x[4*i+j]**2 for j in range(4)) for i in range(3)]; S=sum(r)
    if not S: raise ValueError('Zero field has no normalized value')
    out=sum(a[i]*r[i]**2/2 for i in range(3))
    for e,(i,j) in enumerate(EDGES):
        re=sum(x[4*i+k]*x[4*j+k]+x[4*i+k+1]*x[4*j+k+1] for k in (0,2)); im=sum(x[4*i+k]*x[4*j+k+1]-x[4*i+k+1]*x[4*j+k] for k in (0,2)); out+=b[e]*r[i]*r[j]+c[e]*(re**2+im**2)+d[e]*(re**2-im**2)
    return out/S**2

def find_upper_witness(p,starts=6,seed=2026101605):
    """Only search is floating point; the returned upper bound is exact."""
    from scipy.optimize import minimize
    poly=field_polynomial(**p); mons=np.array(list(poly)); coeff=np.array([float(v) for v in poly.values()]); rng=np.random.default_rng(seed); flags=[]; best=None
    def func(x):
        t=x[mons]; terms=coeff*np.prod(t,axis=1); val=terms.sum(); g=np.zeros(12)
        for k in range(4): np.add.at(g,mons[:,k],coeff*np.prod(np.delete(t,k,axis=1),axis=1))
        S=x@x; return val/S**2,g/S**2-4*val*x/S**3
    for _ in range(starts):
        out=minimize(func,rng.normal(size=12),jac=True,method='BFGS',options={'gtol':1e-10,'maxiter':1200}); flags.append(dict(success=bool(out.success),message=str(out.message),value=float(out.fun)))
        if best is None or out.fun<best.fun: best=out
    y=best.x/np.linalg.norm(best.x); vals=[F(int(round(t*10**9)),10**9) for t in y]; upper=exact_field_value(p,vals)
    return dict(coordinates=list(map(str,vals)),upper_bound=str(upper),starts=starts,optimizer_flags=flags)
