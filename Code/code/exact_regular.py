"""Exact rational implementation of the proposed regular-case full-field reduction.

This checkpoint FAILS CLOSED on singular systems. It is not yet an all-coupling
certified implementation. No neutral-dominance/sign restriction is imposed.
Edge order is (12,23,31). Binary floats are deliberately rejected.
"""
from fractions import Fraction as F
from itertools import combinations, product
from typing import Sequence

EDGES=((0,1),(1,2),(2,0))
ZERO=F(0)

def rational(x):
    if isinstance(x,float): raise TypeError('Use an integer, Fraction or exact decimal/rational string, not binary float')
    return F(x)

def dot(x,y): return sum((a*b for a,b in zip(x,y)),ZERO)
def matvec(A,x): return [dot(row,x) for row in A]
def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]

def solve(A,b):
    n=len(A); M=[[rational(x) for x in row]+[rational(bi)] for row,bi in zip(A,b)]
    if any(len(row)!=n+1 for row in M): raise ValueError('Expected a square matrix')
    for k in range(n):
        p=next((i for i in range(k,n) if M[i][k]),None)
        if p is None: raise ValueError('Singular linear system')
        M[k],M[p]=M[p],M[k]; z=M[k][k]; M[k]=[x/z for x in M[k]]
        for i in range(n):
            if i!=k:
                z=M[i][k]
                if z: M[i]=[x-z*y for x,y in zip(M[i],M[k])]
    return [row[-1] for row in M]

def inverse(A):
    n=len(A); columns=[solve(A,[F(i==j) for i in range(n)]) for j in range(n)]
    return [list(row) for row in zip(*columns)]

def branches(diag,off,delta):
    for k in (1,2,3):
        for I in combinations(range(3),k):
            for tail in product((1,-1),repeat=k-1):
                signs=[1,1,1]
                for i,e in zip(I,(1,*tail)): signs[i]=e
                Q=[[ZERO]*3 for _ in range(3)]
                for i in range(3): Q[i][i]=diag[i]
                for e,(i,j) in enumerate(EDGES): Q[i][j]=Q[j][i]=off[e]+delta[e]*signs[i]*signs[j]
                yield I,Q,'collinear',signs
    if all(delta):
        Q=[[ZERO]*3 for _ in range(3)]
        h=[delta[0]*delta[2]/delta[1],delta[0]*delta[1]/delta[2],delta[1]*delta[2]/delta[0]]
        for i in range(3): Q[i][i]=diag[i]-h[i]
        for e,(i,j) in enumerate(EDGES): Q[i][j]=Q[j][i]=off[e]
        yield (0,1,2),Q,'triangle',None

def phase_cosines(s,kind,signs,delta):
    if kind=='collinear': return [F(signs[i]*signs[j]) for i,j in EDGES]
    g=[abs(delta[0]*delta[2]),abs(delta[0]*delta[1]),abs(delta[1]*delta[2])]; q=[g[i]*s[i] for i in range(3)]
    if any(x<=0 for x in q) or 2*max(q)>sum(q): return None
    sign=lambda x: 1 if x>0 else -1
    sigma=sign(delta[0]*delta[1]*delta[2]); eps=[1,sigma*sign(delta[0]),sigma*sign(delta[2])]
    return [eps[i]*eps[j]*(q[3-i-j]**2-q[i]**2-q[j]**2)/(2*q[i]*q[j]) for i,j in EDGES]

def bordered(Q,I,w):
    return [[Q[i][j] for j in I]+[w[i]] for i in I]+[[F(1)]*len(I)+[ZERO]]

def enumerate_regular(a:Sequence,b:Sequence,c:Sequence,d:Sequence):
    a,b,c,d=([rational(x) for x in v] for v in (a,b,c,d))
    if any(len(v)!=3 for v in (a,b,c,d)): raise ValueError('Each coefficient array must have length three')
    candidates=[]; unresolved=[]; systems=0
    # Neutral rank-one branch; includes all phase-stationary types, not just minima.
    nd=[x/2 for x in d]
    for I,Q,kind,signs in branches([x/2 for x in a],[(b[i]+c[i])/2 for i in range(3)],nd):
        systems+=1
        try: z=solve(bordered(Q,I,[-F(1)]*3),[ZERO]*len(I)+[F(1)])
        except ValueError: unresolved.append(('neutral',I,kind,signs)); continue
        s=[ZERO]*3
        for i,x in zip(I,z): s[i]=x
        if any(s[i]<=0 for i in I): continue
        cos=phase_cosines(s,kind,signs,nd)
        if cos is None: continue
        candidates.append(dict(sector='neutral',support=I,branch=kind,signs=signs,norms=s,cosines=cos,value=z[-1]))
    A=[[ZERO]*3 for _ in range(3)]
    for i in range(3): A[i][i]=a[i]
    for e,(i,j) in enumerate(EDGES): A[i][j]=A[j][i]=b[e]
    alpha=[c[i]+d[i] for i in range(3)]; beta=[c[i]-d[i] for i in range(3)]
    try:
        if not all(alpha) or not all(beta): raise ValueError('Singular overlap block')
        R=inverse(A)
    except ValueError:
        return dict(status='incomplete_singular',unresolved=unresolved+[('charged','singular_L')],candidates=candidates,systems=systems,minimum=None)
    w=matvec(R,[F(1)]*3); traceD=sum(w)
    if not all(w) and traceD:
        r=[x/traceD for x in w]
        if min(r)>=0: candidates.append(dict(sector='charged_nu_zero',support=(),branch='diagonal',norms=r,value=1/(2*traceD)))
    delta=[(1/alpha[i]-1/beta[i])/2 for i in range(3)]
    off=[R[i][j]+(1/alpha[e]+1/beta[e])/2 for e,(i,j) in enumerate(EDGES)]
    for I,Q,kind,signs in branches([R[i][i] for i in range(3)],off,delta):
        systems+=1
        try: z=solve(bordered(Q,I,w),[ZERO]*len(I)+[F(1)])
        except ValueError: unresolved.append(('charged',I,kind,signs)); continue
        s=[ZERO]*3
        for i,x in zip(I,z): s[i]=x
        u=z[-1]
        if any(s[i]<=0 for i in I): continue
        cos=phase_cosines(s,kind,signs,delta)
        if cos is None: continue
        r=matvec(R,s); T=dot(w,s)+u*traceD
        if not T:
            # Trace of E+uD vanishes, so this branch cannot yield trace-one K.
            continue
        trE2=dot(r,r)+sum((s[i]*s[j]*((1/alpha[e]**2+1/beta[e]**2)+(1/alpha[e]**2-1/beta[e]**2)*cos[e]) for e,(i,j) in enumerate(EDGES)),ZERO)
        trK2=(trE2+2*u*dot(r,w)+u*u*dot(w,w))/(T*T)
        if trK2>1: continue
        candidates.append(dict(sector='charged',support=I,branch=kind,signs=signs,kernel_norms=s,kernel_cosines=cos,norms=[(r[i]+u*w[i])/T for i in range(3)],u=u,T=T,trK2=trK2,value=u/(2*T)))
    best=min(candidates,key=lambda x:x['value']) if candidates else None
    return dict(status='incomplete_singular' if unresolved else 'regular_case_exact_minimum',unresolved=unresolved,candidates=candidates,systems=systems,minimum=best if not unresolved else None)

if __name__=='__main__':
    import json
    p=dict(a=['31/5','36/5','41/5'],b=['-19/5','-37/10','-39/10'],c=['1','6/5','4/5'],d=['1/10','3/25','3/20'])
    result=enumerate_regular(**p); print(json.dumps(result,default=str,indent=2))
