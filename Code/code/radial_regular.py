"""Second regular-domain formulation: solve in original doublet norms.

For full kernel support, Q_dual=A^{-1}+M implies
(A+M^{-1})r=mu*1, sum r=1. The code never inverts A to obtain candidates.
The 3 orthogonal-pair and 6 parallel-pair/orthogonal-singleton candidates
are built directly from physical angular configurations. Neutral candidates
are recomputed using independent SymPy linear algebra. All singular systems
fail closed; full_infinitesimal.py remains the all-coupling implementation.
"""
from fractions import Fraction as F
from itertools import combinations,product
import sympy as sp
from exact_witness import Quad,phase_units,verify_certificate,determinant
EDGES=((0,1),(1,2),(2,0))

def solve_radial(B,I):
    n=len(I); mat=B.extract(I,I).row_join(-sp.ones(n,1)).col_join(sp.ones(1,n).row_join(sp.zeros(1,1)))
    if mat.det()==0: raise ArithmeticError('Singular radial stationary system')
    x=mat.inv()*sp.Matrix([0]*n+[1]); r=[F(0)]*3
    for i,v in zip(I,x[:-1]): r[i]=F(v)
    return r,F(x[-1])

def all_phase_data(diag,off,delta):
    for size in (1,2,3):
        for I in combinations(range(3),size):
            for tail in product((1,-1),repeat=size-1):
                signs=[1]*3
                for i,e in zip(I,(1,*tail)): signs[i]=e
                Q=sp.diag(*diag)
                for e,(i,j) in enumerate(EDGES): Q[i,j]=Q[j,i]=off[e]+delta[e]*signs[i]*signs[j]
                yield I,Q,'collinear',signs
    if not all(delta): raise ArithmeticError('Non-generic phase graph')
    Q=sp.diag(*diag); h=[delta[0]*delta[2]/delta[1],delta[0]*delta[1]/delta[2],delta[1]*delta[2]/delta[0]]
    for i in range(3): Q[i,i]-=h[i]
    for e,(i,j) in enumerate(EDGES): Q[i,j]=Q[j,i]=off[e]
    yield (0,1,2),Q,'triangle',None

def radial_minimum(a,b,c,d):
    if any(len(v)!=3 for v in (a,b,c,d)): raise ValueError('Each coefficient array must have length three')
    if any(isinstance(x,float) for v in (a,b,c,d) for x in v): raise TypeError('Binary floats are not accepted')
    p={k:[F(v) for v in vals] for k,vals in dict(a=a,b=b,c=c,d=d).items()}; a,b,c,d=[p[k] for k in ('a','b','c','d')]; A=sp.diag(*a); candidates=[]; systems=0; unresolved=[]
    for e,(i,j) in enumerate(EDGES): A[i,j]=A[j,i]=b[e]
    if A.det()==0: return dict(status='incomplete_singular',systems=0,minimum=None,unresolved=['Singular radial matrix A'],candidates=[])
    def store(kind,r,mu,**extra): candidates.append(dict(kind=kind,radii=list(map(str,r)),value=mu/2,**extra))
    try:
        # Neutral (14 branches), with direct phase realizability.
        for I,Q,kind,signs in all_phase_data([x/2 for x in a],[(b[i]+c[i])/2 for i in range(3)],[x/2 for x in d]):
            systems+=1
            try: r,v=solve_radial(Q,I)
            except ArithmeticError as e: unresolved.append(str(e)); continue
            if any(r[i]<=0 for i in I): continue
            try: phase_units(r,kind,signs,[x/2 for x in d])
            except ValueError: continue
            store('neutral_'+kind,r,2*v,support=I)
        # Three physical orthogonal pairs (one doublet absent).
        for I in combinations(range(3),2):
            systems+=1
            try: r,mu=solve_radial(A,I)
            except ArithmeticError as e: unresolved.append(str(e)); continue
            if all(r[i]>0 for i in I): store('charged_orthogonal_pair',r,mu,support=I)
        # Six parallel-pair / orthogonal-singleton branches.
        for e,(i,j) in enumerate(EDGES):
            for phase_sign in (1,-1):
                B=A.copy(); B[i,j]+=c[e]+phase_sign*d[e]; B[j,i]=B[i,j]; systems+=1
                try: r,mu=solve_radial(B,(0,1,2))
                except ArithmeticError as ex: unresolved.append(str(ex)); continue
                if min(r)>0: store('charged_split',r,mu,pair=(i,j),phase_sign=phase_sign)
        # Five full-kernel branches: radial effective matrix A + inverse(M).
        alpha=[c[i]+d[i] for i in range(3)]; beta=[c[i]-d[i] for i in range(3)]
        if not all(alpha) or not all(beta): raise ArithmeticError('Singular overlap block')
        delta=[(1/alpha[i]-1/beta[i])/2 for i in range(3)]; kappa=[(1/alpha[i]+1/beta[i])/2 for i in range(3)]
        for I,M,kind,signs in all_phase_data([F(0)]*3,kappa,delta):
            if len(I)!=3: continue
            systems+=1
            if M.det()==0: unresolved.append('Singular angular matrix M'); continue
            Minv=M.inv(); B=A+Minv
            try: r,mu=solve_radial(B,I)
            except ArithmeticError as ex: unresolved.append(str(ex)); continue
            if min(r)<=0: continue
            zeta=[F(v) for v in Minv*sp.Matrix(r)]; Z=sum(zeta)
            if not Z: continue
            s=[v/Z for v in zeta]
            if min(s)<=0: continue
            try: h,z=phase_units(s,kind,signs,delta)
            except ValueError: continue
            H=[[Quad(0,0,h) for _ in range(3)] for _ in range(3)]
            for i in range(3): H[i][i]=Quad(s[i]*r[i],0,h)
            for e,(i,j) in enumerate(EDGES): H[i][j]=-Z*s[i]*s[j]*(kappa[e]+delta[e]*z[i].conjugate()*z[j]); H[j][i]=H[i][j].conjugate()
            minors=[H[i][i].a for i in range(3)]+[H[i][i].a*H[j][j].a-H[i][j].norm2() for i,j in EDGES]
            if min(minors)<0: continue
            if not determinant(H).zero(): raise ArithmeticError('Direct construction lost its null vector')
            cert=dict(coefficients={k:list(map(str,v))for k,v in p.items()},sector='charged',branch=kind,s=list(map(str,s)),radicand=str(h),phase_units=[x.pair() for x in z],H=[[x.pair() for x in row]for row in H],value=str(mu/2),u=str(-mu/Z),T=str(-1/Z)); checked=verify_certificate(cert)
            store('charged_full_'+kind,r,mu,kernel_norms=list(map(str,s)),certificate=cert,certificate_verified=checked)
    except ArithmeticError as e: unresolved.append(str(e))
    best=min(candidates,key=lambda r:r['value']) if candidates else None
    return dict(status='incomplete_singular' if unresolved else 'regular_radial_result',systems=systems,minimum=best if not unresolved else None,unresolved=unresolved,candidates=candidates)
