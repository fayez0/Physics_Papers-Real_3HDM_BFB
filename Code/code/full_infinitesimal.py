"""Experimental all-coupling extension using a positive infinitesimal.

Proposed theorem and implementation require an independent audit. This is NOT
an acceptance-ready publication claim. Arithmetic is exact over Q(eps); the
order is the sign for eps -> 0+, never a small floating-point substitution.

Potential/edge convention: same as exact_regular.py. All rational inputs,
including singular and quartic-flat inputs, are accepted. Binary floats fail.
"""
from fractions import Fraction
from itertools import combinations, product
from sympy.polys.fields import field
from sympy.polys.domains import QQ

RF, EPS = field('eps', QQ)
Z,O=RF.zero,RF.one
EDGES=((0,1),(1,2),(2,0))
PERTURBATION={'a':[11,6,19],'b':[3,-2,3],'c':[6,11,3],'d':[1,1,2]}

def scalar(x):
    if isinstance(x,float): raise TypeError('Binary floats are not accepted; supply rational or exact decimal strings')
    r=Fraction(x); return RF(QQ(r.numerator,r.denominator))

def leading_at_zero(x):
    if not x: return None
    n=min(x.numer.to_dict().items()); d=min(x.denom.to_dict().items())
    return n[0][0]-d[0][0],n[1]/d[1]

def sgn(x):
    if not x: return 0
    a=leading_at_zero(x)[1]; return 1 if a>0 else -1

def absolute(x): return x if sgn(x)>=0 else -x

def limit_zero(x):
    if not x: return Fraction(0)
    power,coefficient=leading_at_zero(x)
    if power<0: raise ArithmeticError('Divergent value on an allegedly physical branch')
    if power>0: return Fraction(0)
    return Fraction(int(coefficient.numerator),int(coefficient.denominator))

def dot(x,y): return sum((a*b for a,b in zip(x,y)),Z)
def mv(A,x): return [dot(row,x) for row in A]

def solve(A,b):
    n=len(A); M=[row[:]+[bi] for row,bi in zip(A,b)]
    for k in range(n):
        pivot=next((j for j in range(k,n) if M[j][k]),None)
        if pivot is None: raise ArithmeticError('Regularizer failed: identically singular system')
        M[k],M[pivot]=M[pivot],M[k]; scale=M[k][k]; M[k]=[x/scale for x in M[k]]
        for j in range(n):
            if j==k: continue
            scale=M[j][k]
            if scale: M[j]=[x-scale*y for x,y in zip(M[j],M[k])]
    return [row[-1] for row in M]

def inverse(A):
    columns=[solve(A,[O if i==j else Z for i in range(3)]) for j in range(3)]
    return [list(row) for row in zip(*columns)]

def bordered(Q,I,w): return [[Q[i][j] for j in I]+[w[i]] for i in I]+[[O]*len(I)+[Z]]

def branches(diag,off,delta):
    for size in (1,2,3):
        for I in combinations(range(3),size):
            for tail in product((1,-1),repeat=size-1):
                signs=[1,1,1]
                for i,v in zip(I,(1,*tail)): signs[i]=v
                Q=[[Z]*3 for _ in range(3)]
                for i in range(3): Q[i][i]=diag[i]
                for e,(i,j) in enumerate(EDGES): Q[i][j]=Q[j][i]=off[e]+delta[e]*signs[i]*signs[j]
                yield I,Q,'collinear',signs
    # All delta entries are nonzero rational functions for the fixed regularizer.
    if not all(delta): raise ArithmeticError('Regularizer failed to make the phase graph generic')
    Q=[[Z]*3 for _ in range(3)]
    h=[delta[0]*delta[2]/delta[1],delta[0]*delta[1]/delta[2],delta[1]*delta[2]/delta[0]]
    for i in range(3): Q[i][i]=diag[i]-h[i]
    for e,(i,j) in enumerate(EDGES): Q[i][j]=Q[j][i]=off[e]
    yield (0,1,2),Q,'triangle',None

def phase_cosines(s,kind,signs,delta):
    if kind=='collinear': return [RF(signs[i]*signs[j]) for i,j in EDGES]
    g=[absolute(delta[0]*delta[2]),absolute(delta[0]*delta[1]),absolute(delta[1]*delta[2])]; q=[g[i]*s[i] for i in range(3)]
    if any(sgn(x)<=0 for x in q) or any(sgn(sum(q,Z)-2*x)<0 for x in q): return None
    sigma=sgn(delta[0]*delta[1]*delta[2]); eta=[1,sigma*sgn(delta[0]),sigma*sgn(delta[2])]
    return [eta[i]*eta[j]*(q[3-i-j]**2-q[i]**2-q[j]**2)/(2*q[i]*q[j]) for i,j in EDGES]

def full_minimum(a,b,c,d,keep_functions=False):
    inputs=dict(a=a,b=b,c=c,d=d)
    if any(len(x)!=3 for x in inputs.values()): raise ValueError('Expected four arrays of length three')
    p={name:[scalar(x)+EPS*PERTURBATION[name][i] for i,x in enumerate(xs)] for name,xs in inputs.items()}
    a,b,c,d=(p[name] for name in ('a','b','c','d')); candidates=[]; systems=0
    def store(sector,I,kind,signs,value,**extras):
        row=dict(sector=sector,support=I,branch=kind,signs=signs,value_limit=limit_zero(value))
        if keep_functions: row.update(value_function=str(value),**{k:[str(x) for x in v] if isinstance(v,list) else str(v) for k,v in extras.items()})
        candidates.append(row)
    nd=[x/2 for x in d]
    for I,Q,kind,signs in branches([x/2 for x in a],[(b[i]+c[i])/2 for i in range(3)],nd):
        systems+=1; z=solve(bordered(Q,I,[-O]*3),[Z]*len(I)+[O]); s=[Z]*3
        for i,x in zip(I,z): s[i]=x
        if any(sgn(s[i])<=0 for i in I): continue
        cos=phase_cosines(s,kind,signs,nd)
        if cos is not None: store('neutral',I,kind,signs,z[-1],norms=s)
    A=[[Z]*3 for _ in range(3)]
    for i in range(3): A[i][i]=a[i]
    for e,(i,j) in enumerate(EDGES): A[i][j]=A[j][i]=b[e]
    R=inverse(A); w=mv(R,[O]*3); d0=sum(w,Z); alpha=[c[i]+d[i] for i in range(3)]; beta=[c[i]-d[i] for i in range(3)]
    if not all(alpha) or not all(beta) or not all(w): raise ArithmeticError('Regularizer failed in inverse map')
    delta=[(1/alpha[i]-1/beta[i])/2 for i in range(3)]
    off=[R[i][j]+(1/alpha[e]+1/beta[e])/2 for e,(i,j) in enumerate(EDGES)]
    for I,Q,kind,signs in branches([R[i][i] for i in range(3)],off,delta):
        systems+=1; z=solve(bordered(Q,I,w),[Z]*len(I)+[O]); s=[Z]*3
        for i,x in zip(I,z): s[i]=x
        u=z[-1]
        if any(sgn(s[i])<=0 for i in I): continue
        cos=phase_cosines(s,kind,signs,delta)
        if cos is None: continue
        T=dot(w,s)+u*d0
        if not T: continue
        r=mv(R,s)
        trE2=dot(r,r)+sum((s[i]*s[j]*((1/alpha[e]**2+1/beta[e]**2)+(1/alpha[e]**2-1/beta[e]**2)*cos[e]) for e,(i,j) in enumerate(EDGES)),Z)
        trK2=(trE2+2*u*dot(r,w)+u*u*dot(w,w))/T**2
        if sgn(O-trK2)<0: continue
        store('charged',I,kind,signs,u/(2*T),kernel_norms=s,trK2=trK2)
    if not candidates: raise ArithmeticError('No physical branch survived; do not classify the input')
    best=min(candidates,key=lambda x:x['value_limit']); m=best['value_limit']; classification='strict_quartic_positive' if m>0 else 'quartic_nonnegative_with_flat_direction' if m==0 else 'quartic_unbounded_below'
    return dict(status='experimental_full_domain_result',classification=classification,minimum=m,best_branch_limit=best,systems=systems,regularizer=PERTURBATION,candidates=candidates)

if __name__=='__main__':
    import json,time,sys
    mode=sys.argv[1] if len(sys.argv)>1 else 'zero'
    examples={'zero':([0]*3,[0]*3,[0]*3,[0]*3),'rank_flat':([1]*3,[-1]*3,[2]*3,[0]*3),'triangle':([3]*3,[-2]*3,[1]*3,['1/2']*3),'positive':([1]*3,[0]*3,[0]*3,[0]*3)}
    t=time.monotonic(); result=full_minimum(*examples[mode]); result['elapsed_seconds']=time.monotonic()-t; print(json.dumps(result,default=str,indent=2))
