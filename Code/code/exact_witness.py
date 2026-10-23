"""Exact physical Gram certificates, independent of the scalar trace-square filter.

A certificate stores H=B^dagger K B over Q(i sqrt(h)), h>=0, with invertible
B_ii=sqrt(s_i)*sqrt(z_i) for active kernel/neutral components and B_ii=1
otherwise. Thus PSD/rank follow from H's principal minors, without rounded
phases, eigenvalues, or the solver's trace-square expression. Physical value
and full stationarity are recomputed directly from the original couplings.
Certificates prove candidate realizability, NOT completeness of the candidate list.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
import sympy as sp

EDGES=((0,1),(1,2),(2,0))

@dataclass(frozen=True)
class Quad:
    a: F
    b: F
    h: F
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a)); object.__setattr__(self,'h',F(self.h)); object.__setattr__(self,'b',F(self.b) if self.h else F(0))
        if self.h<0: raise ValueError('Expected real nonnegative radicand')
    def cast(self,o):
        if isinstance(o,Quad):
            if o.h!=self.h: raise ValueError('Quadratic fields differ')
            return o
        return Quad(F(o),F(0),self.h)
    def __add__(self,o):
        o=self.cast(o); return Quad(self.a+o.a,self.b+o.b,self.h)
    __radd__=__add__
    def __neg__(self): return Quad(-self.a,-self.b,self.h)
    def __sub__(self,o): return self+-self.cast(o)
    def __rsub__(self,o): return self.cast(o)+-self
    def __mul__(self,o):
        o=self.cast(o); return Quad(self.a*o.a-self.h*self.b*o.b,self.a*o.b+self.b*o.a,self.h)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=self.cast(o); den=o.a*o.a+self.h*o.b*o.b
        if not den: raise ZeroDivisionError
        return self*Quad(o.a/den,-o.b/den,self.h)
    def conjugate(self): return Quad(self.a,-self.b,self.h)
    def norm2(self): return self.a*self.a+self.h*self.b*self.b
    def pair(self): return [str(self.a),str(self.b)]
    def zero(self): return self.a==0 and self.b==0

def phase_units(s,branch,signs,delta):
    if branch=='collinear': return F(0),[Quad(F(v),0,0) for v in signs]
    g=[abs(delta[0]*delta[2]),abs(delta[0]*delta[1]),abs(delta[1]*delta[2])]; q=[g[i]*s[i] for i in range(3)]
    if min(q)<=0 or 2*max(q)>sum(q): raise ValueError('Unphysical phase triangle')
    numerator=q[2]**2-q[0]**2-q[1]**2; h=4*q[0]**2*q[1]**2-numerator**2
    z=[Quad(1,0,h),Quad(numerator/(2*q[0]*q[1]),1/(2*q[0]*q[1]),h)]
    z.append(-(q[0]+q[1]*z[1])/q[2]); sign=lambda x:1 if x>0 else -1; sigma=sign(delta[0]*delta[1]*delta[2]); eta=[1,sigma*sign(delta[0]),sigma*sign(delta[2])]
    return h,[eta[i]*z[i] for i in range(3)]

def determinant(H):
    out=H[0][0].cast(0)
    for perm in permutations(range(3)):
        inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)); term=H[0][0].cast((-1)**inversions)
        for i in range(3): term=term*H[i][perm[i]]
        out=out+term
    return out

def construct_certificate(p,candidate):
    p={k:[F(x) for x in v] for k,v in p.items()}; a,b,c,d=[p[k] for k in ('a','b','c','d')]; sector=candidate['sector']; s=[F(x) for x in candidate['norms' if sector=='neutral' else 'kernel_norms']]
    if sector=='neutral': delta=[x/2 for x in d]
    elif sector=='charged': delta=[(1/(c[i]+d[i])-1/(c[i]-d[i]))/2 for i in range(3)]
    else: raise ValueError('Use generic rank-one or rank-two candidate')
    h,z=phase_units(s,candidate['branch'],candidate['signs'],delta); Z=Quad(0,0,h); H=[[Z for _ in range(3)] for _ in range(3)]; t=[x if x else F(1) for x in s]
    if sector=='neutral':
        H=[[Quad(s[i]*s[j],0,h) for j in range(3)] for i in range(3)]
    else:
        A=sp.diag(*a)
        for e,(i,j) in enumerate(EDGES): A[i,j]=A[j,i]=b[e]
        R=A.inv(); w=[F(x) for x in R*sp.ones(3,1)]; rE=[F(x) for x in R*sp.Matrix(s)]; u=F(candidate['u']); T=F(candidate['T'])
        for i in range(3): H[i][i]=Quad(t[i]*(rE[i]+u*w[i])/T,0,h)
        for e,(i,j) in enumerate(EDGES):
            kappa=(1/(c[e]+d[e])+1/(c[e]-d[e]))/2; phase=z[i].conjugate()*z[j]; H[i][j]=(s[i]*s[j]/T)*(kappa+delta[e]*phase); H[j][i]=H[i][j].conjugate()
    return dict(coefficients={k:list(map(str,v)) for k,v in p.items()},sector=sector,branch=candidate['branch'],s=list(map(str,s)),radicand=str(h),phase_units=[x.pair() for x in z],H=[[x.pair() for x in row] for row in H],value=str(candidate['value']),u=str(candidate.get('u',0)),T=str(candidate.get('T',1)))

def verify_certificate(cert):
    def require(condition, message):
        if not condition: raise ValueError('Invalid Gram certificate: '+message)
    require(cert.get('sector') in ('neutral','charged'), 'unknown sector')
    require(len(cert['s'])==3 and len(cert['phase_units'])==3 and len(cert['H'])==3 and all(len(row)==3 for row in cert['H']), 'wrong dimensions')
    p={k:[F(x) for x in v] for k,v in cert['coefficients'].items()}; a,b,c,d=[p[k] for k in ('a','b','c','d')]; s=[F(x) for x in cert['s']]; h=F(cert['radicand']); z=[Quad(*pair,h) for pair in cert['phase_units']]; H=[[Quad(*pair,h) for pair in row] for row in cert['H']]; t=[x if x else F(1) for x in s]; target=F(cert['value']); active=[int(x>0) for x in s]; Z=Quad(0,0,h)
    require(sum(s)==1 and min(s)>=0 and all(zi.norm2()==1 for zi in z), 'kernel normalization or phase unit norm')
    require(all((H[i][j]-H[j][i].conjugate()).zero() for i in range(3) for j in range(3)), 'non-Hermitian matrix')
    r=[H[i][i].a/t[i] for i in range(3)]; require(sum(r)==1, 'Gram trace not one')
    minors=[H[i][i].a for i in range(3)]+[H[i][i].a*H[j][j].a-H[i][j].norm2() for i,j in EDGES]; det=determinant(H)
    require(min(minors)>=0 and det.zero(), 'PSD or rank constraint failed')
    v=sum(a[i]*r[i]*r[i]/2 for i in range(3)); Lh=[[Z for _ in range(3)] for _ in range(3)]
    A=[[F(0)]*3 for _ in range(3)]
    for i in range(3): A[i][i]=a[i]
    for e,(i,j) in enumerate(EDGES):
        A[i][j]=A[j][i]=b[e]; phase=z[i].conjugate()*z[j]; reK2=(H[i][j]*H[i][j]*phase.conjugate()).a/(t[i]*t[j]); absK2=H[i][j].norm2()/(t[i]*t[j]); v+=b[e]*r[i]*r[j]+c[e]*absK2+d[e]*reK2
        Lh[i][j]=c[e]*H[i][j]+d[e]*H[i][j].conjugate()*phase; Lh[j][i]=Lh[i][j].conjugate()
    for i in range(3): Lh[i][i]=Quad(t[i]*sum(A[i][j]*r[j] for j in range(3)),0,h)
    require(v==target, 'original potential value mismatch')
    if cert['sector']=='charged':
        require(all(sum((H[i][j]*active[j] for j in range(3)),Z).zero() for i in range(3)), 'charged null vector failed'); T=F(cert['T']); mu=2*target
        require(all((Lh[i][j]-(mu*t[i] if i==j else 0)-s[i]*s[j]/T).zero() for i in range(3) for j in range(3)), 'charged stationarity failed')
    else:
        require(all((sum((Lh[i][j]*active[j] for j in range(3)),Z)-2*target*s[i]).zero() for i in range(3)), 'neutral stationarity failed')
    return dict(physical=True,stationarity_exact=True,normalized_value=str(v),principal_minors=list(map(str,minors)),determinant=det.pair(),radii=list(map(str,r)),extension_radicand=str(h))
