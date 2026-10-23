"""All-coupling original-norm formulation; no inverse of A or inverse-gradient solver.

Every input is deformed by the fixed, exactly checked auxiliary coefficients.
Branches are physical orthogonal/split configurations and five full-kernel
angular forms. Small stationary systems use explicit adjugates, not an inverse
of the radial coefficient matrix. This is a second implementation of the
research theorem, not a stand-alone proof certificate of nonnegativity.
"""
from fractions import Fraction as F
from itertools import combinations,product
from sympy.polys.fields import field
from sympy.polys.domains import QQ
RF,E=field('eps',QQ); Z,O=RF.zero,RF.one
EDGES=((0,1),(1,2),(2,0)); P0={'a':[11,6,19],'b':[3,-2,3],'c':[6,11,3],'d':[1,1,2]}

def q(x):
    if isinstance(x,float):raise TypeError('Binary floats are not accepted')
    x=F(x);return RF(QQ(x.numerator,x.denominator))

def sign(x):
    if not x:return 0
    a=min(x.numer.terms(),key=lambda t:t[0])[1];b=min(x.denom.terms(),key=lambda t:t[0])[1];return 1 if a*b>0 else -1

def absq(x):return x if sign(x)>=0 else -x

def limit(x):
    if not x:return F(0)
    (i,),a=min(x.numer.terms(),key=lambda t:t[0]);(j,),b=min(x.denom.terms(),key=lambda t:t[0])
    if i<j:raise ArithmeticError('Divergent physical value')
    return F(int(a.numerator),int(a.denominator))/F(int(b.numerator),int(b.denominator)) if i==j else F(0)

def matrix(diag,off):
    A=[[Z for _ in range(3)] for _ in range(3)]
    for i in range(3):A[i][i]=diag[i]
    for e,(i,j) in enumerate(EDGES):A[i][j]=A[j][i]=off[e]
    return A

def adjugate(B):
    a,b,c=B[0];_,d,e=B[1];_,_,f=B[2]
    return [[d*f-e*e,c*e-b*f,b*e-c*d],[c*e-b*f,a*f-c*c,b*c-a*e],[b*e-c*d,b*c-a*e,a*d-b*b]]

def det3(B):
    C=adjugate(B);return sum((B[0][j]*C[j][0] for j in range(3)),Z)

def stationary(B,I):
    r=[Z,Z,Z]
    if len(I)==1:r[I[0]]=O;return r,B[I[0]][I[0]],-O
    if len(I)==2:
        i,j=I;D=B[i][i]+B[j][j]-2*B[i][j]
        if not D:raise ArithmeticError('Identically singular radial bordered system')
        r[i]=(B[j][j]-B[i][j])/D;r[j]=(B[i][i]-B[i][j])/D;return r,(B[i][i]*B[j][j]-B[i][j]**2)/D,-D
    C=adjugate(B);w=[sum(row,Z) for row in C];D=sum(w,Z)
    if not D:raise ArithmeticError('Identically singular full radial bordered system')
    return [x/D for x in w],det3(B)/D,-D

def phases(s,delta,kind,signs=None):
    if kind=='collinear':return [RF(signs[i]*signs[j]) for i,j in EDGES]
    g=[absq(delta[0]*delta[2]),absq(delta[0]*delta[1]),absq(delta[1]*delta[2])];Q=[g[i]*s[i] for i in range(3)];S=sum(Q,Z)
    if any(sign(x)<=0 or sign(S-2*x)<0 for x in Q):return None
    sg=sign(delta[0]*delta[1]*delta[2]);eta=[1,sg*sign(delta[0]),sg*sign(delta[2])]
    return [eta[i]*eta[j]*(Q[3-i-j]**2-Q[i]**2-Q[j]**2)/(2*Q[i]*Q[j]) for i,j in EDGES]

def angular_inverse(c,d,kind,signs=None):
    if kind=='collinear':
        gamma=[c[e]+signs[i]*signs[j]*d[e] for e,(i,j) in enumerate(EDGES)]
        diag=[-gamma[0]*gamma[2]/(2*gamma[1]),-gamma[0]*gamma[1]/(2*gamma[2]),-gamma[1]*gamma[2]/(2*gamma[0])]
        return matrix(diag,[v/2 for v in gamma])
    f=[c[e]**2-d[e]**2 for e in range(3)];D=d[0]*d[1]*d[2];H=D*D+2*c[0]*c[1]*c[2]*D-sum((c[e]**2*d[(e+1)%3]**2*d[(e+2)%3]**2 for e in range(3)),Z)
    if not H:raise ArithmeticError('Identically singular angular triangle matrix')
    diag=[-f[0]*f[2]*D/H,-f[0]*f[1]*D/H,-f[1]*f[2]*D/H];off=[]
    for e in range(3):
        j,k=[z for z in range(3) if z!=e];off.append(f[e]*d[j]*d[k]*(c[j]*c[k]*d[e]-c[e]*d[j]*d[k])/H)
    return matrix(diag,off)

def run_radial(a,b,c,d,deform=True,retain_functions=False):
    inputs=dict(a=a,b=b,c=c,d=d)
    if any(len(v)!=3 for v in inputs.values()):raise ValueError('Expected four length-three arrays')
    p={k:[q(v)+(E*q(P0[k][i]) if deform else Z) for i,v in enumerate(vals)] for k,vals in inputs.items()};a,b,c,d=[p[k] for k in ('a','b','c','d')];A=matrix(a,b);records=[];candidates=[]
    if any(not x for x in d+[c[e]+d[e] for e in range(3)]+[c[e]-d[e] for e in range(3)]):raise ArithmeticError('Nonregular angular coefficients; use exact deformation')
    def store(label,B,I,checker=None,**meta):
        r,mu,bdet=stationary(B,I);row=dict(label=label,support=list(I),**meta);row['bordered_determinant_at_p0']=str(bdet) if not deform else None
        if retain_functions:row.update(radii=[str(x) for x in r],value_function=str(mu/2),value_degrees=[(mu/2).numer.degree(),(mu/2).denom.degree()])
        if any(sign(r[i])<=0 for i in I):row['rejection']='nonpositive_radii';records.append(row);return
        if checker:
            reason=checker(r)
            if reason:row['rejection']=reason;records.append(row);return
        row.update(accepted=True,minimum_limit=limit(mu/2));records.append(row);candidates.append(row)
    # Rank-one configurations: thirteen phase-collinear supports plus one triangle.
    for size in (1,2,3):
        for I in combinations(range(3),size):
            for tail in product((1,-1),repeat=size-1):
                eps=[1]*3
                for i,t in zip(I,(1,*tail)):eps[i]=t
                B=matrix(a,[b[e]+c[e]+d[e]*eps[i]*eps[j] for e,(i,j) in enumerate(EDGES)]);store('neutral_collinear',B,I,phase_signs=eps)
    B=matrix([a[0]-d[0]*d[2]/d[1],a[1]-d[0]*d[1]/d[2],a[2]-d[1]*d[2]/d[0]],[b[e]+c[e] for e in range(3)])
    store('neutral_triangle',B,(0,1,2),lambda r:'unrealizable_triangle' if phases(r,d,'triangle') is None else None)
    # Three orthogonal pairs and six parallel-pair / orthogonal-singleton branches.
    for I in combinations(range(3),2):store('charged_orthogonal_pair',A,I)
    for e,(i,j) in enumerate(EDGES):
        for t in (1,-1):
            B=[row[:] for row in A];B[i][j]+=c[e]+t*d[e];B[j][i]=B[i][j];store('charged_split',B,(0,1,2),pair=[i,j],phase_sign=t)
    alpha=[c[e]+d[e] for e in range(3)];beta=[c[e]-d[e] for e in range(3)];delta=[(1/alpha[e]-1/beta[e])/2 for e in range(3)]
    # Four full-support collinear kernels and one full-support triangle kernel.
    for kind,eps in [('collinear',[1,*tail]) for tail in product((1,-1),repeat=2)]+[('triangle',None)]:
        G=angular_inverse(c,d,kind,eps);B=[[A[i][j]+G[i][j] for j in range(3)] for i in range(3)]
        def check(r,G=G,kind=kind,eps=eps):
            v=[sum((G[i][j]*r[j] for j in range(3)),Z) for i in range(3)];total=sum(v,Z)
            if not total:return 'zero_kernel_normalization'
            s=[t/total for t in v]
            if any(sign(x)<=0 for x in s):return 'nonpositive_kernel_support'
            cos=phases(s,delta,kind,eps)
            if cos is None:return 'unrealizable_kernel_triangle'
            tr=sum((x*x for x in r),Z)+total**2*sum((s[i]*s[j]*(1/alpha[e]**2+1/beta[e]**2+(1/alpha[e]**2-1/beta[e]**2)*cos[e]) for e,(i,j) in enumerate(EDGES)),Z)
            return 'indefinite_Gram' if sign(O-tr)<0 else None
        store('charged_full_'+kind,B,(0,1,2),check,phase_signs=eps)
    if len(records)!=28 or not candidates:raise ArithmeticError('Incomplete branch list')
    best=min(candidates,key=lambda row:row['minimum_limit']);m=best['minimum_limit'];return dict(status='full_radial_exact_result' if deform else 'direct_regular_result',minimum=m,best=best,classification='strict_quartic_positive' if m>0 else 'quartic_unbounded_below' if m<0 else 'quartic_nonnegative_with_flat_direction',systems=28,records=records)

if __name__=='__main__':
    import json,time
    t=time.monotonic();out=run_radial([3]*3,[-2]*3,[1]*3,['1/2']*3,retain_functions=True);out['elapsed_seconds']=time.monotonic()-t;print(json.dumps(out,default=str,indent=2))
