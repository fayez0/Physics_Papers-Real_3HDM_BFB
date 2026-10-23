"""Standalone exact full-decision transcript checker. Imports NO BFB solver.

Trust boundary: the mathematical 28-branch completeness/deformation theorem,
Python, and SymPy rational-polynomial arithmetic. No minimization, no equation
solver, no floating epsilon, and no eval/exec of certificate text is used.
"""
from fractions import Fraction as F
from itertools import combinations,product
from sympy.polys.fields import field
from sympy.polys.domains import QQ
FIELD,e=field('t',QQ);ZERO,ONE=FIELD.zero,FIELD.one
EDGES=((0,1),(1,2),(2,0))
P0=dict(a=[11,6,19],b=[3,-2,3],c=[6,11,3],d=[1,1,2])

def require(test,message):
    if not test:raise ValueError('Invalid decision certificate: '+message)

def rational(v):
    require(not isinstance(v,(float,bool)),'nonexact input')
    x=F(v);return FIELD(QQ(x.numerator,x.denominator))

def read_function(obj):
    require(isinstance(obj,dict) and set(obj)=={'numerator','denominator'},'malformed rational function')
    def polynomial(values):
        require(isinstance(values,list) and 1<=len(values)<=501,'invalid polynomial size')
        out=ZERO
        for v in reversed(values):out=out*e+rational(v)
        return out
    n,d=polynomial(obj['numerator']),polynomial(obj['denominator']);require(bool(d),'zero polynomial denominator')
    return n/d

def lead(v):
    if not v:return None
    nt=v.numer.to_dict();dt=v.denom.to_dict();ni=min(k[0] for k in nt);di=min(k[0] for k in dt)
    return ni-di,nt[(ni,)]/dt[(di,)]

def sign(v):
    l=lead(v);return 0 if l is None else (1 if l[1]>0 else -1)

def limit(v):
    l=lead(v)
    if l is None:return F(0)
    require(l[0]>=0,'divergent physical value')
    return F(0) if l[0]>0 else F(int(l[1].numerator),int(l[1].denominator))

def dot(v,w):return sum((x*y for x,y in zip(v,w)),ZERO)

def descriptors():
    out=[]
    for sector in ['neutral','charged']:
        for size in [1,2,3]:
            for I in combinations(range(3),size):
                for tail in product([1,-1],repeat=size-1):
                    signs=[1,1,1]
                    for i,sgn in zip(I,(1,*tail)):signs[i]=sgn
                    out.append((sector,list(I),'collinear',signs))
        out.append((sector,[0,1,2],'triangle',None))
    return out

def verify_decision(cert):
    require(cert.get('schema')=='weinberg-full-decision-v1','wrong schema')
    require(cert.get('regularizer')==P0,'uncertified regularizer')
    raw=cert['coefficients'];require(set(raw)==set(P0) and all(len(v)==3 for v in raw.values()),'coefficient dimensions')
    p={k:[rational(v)+e*P0[k][i] for i,v in enumerate(raw[k])] for k in P0};a,b,c,d=[p[k] for k in ['a','b','c','d']]
    # Explicit adjugate: no call to a matrix inverse or equation solver.
    a0,a1,a2=a;b0,b1,b2=b
    det=a0*a1*a2+2*b0*b1*b2-a0*b1*b1-a1*b2*b2-a2*b0*b0
    require(bool(det),'radial map identically singular')
    R=[[ (a1*a2-b1*b1)/det,(b1*b2-b0*a2)/det,(b0*b1-a1*b2)/det],[(b1*b2-b0*a2)/det,(a0*a2-b2*b2)/det,(b0*b2-a0*b1)/det],[(b0*b1-a1*b2)/det,(b0*b2-a0*b1)/det,(a0*a1-b0*b0)/det]]
    w=[sum(row,ZERO) for row in R];alpha=[c[i]+d[i] for i in range(3)];beta=[c[i]-d[i] for i in range(3)]
    require(all(w) and all(alpha) and all(beta) and all(d),'identically singular map/phase coefficients')
    rows=cert['rows'];expected=descriptors();require(len(rows)==len(expected)==28,'all 28 branches must be present')
    values=[];counts={};rejects={};details=[]
    for row,(sector,I,kind,eps) in zip(rows,expected):
        require((row.get('sector'),row.get('support'),row.get('kind'),row.get('signs'))==(sector,I,kind,eps),'branch missing, duplicated, reordered or altered')
        require(len(row['s'])==3,'wrong squared-component dimension');s=[read_function(v) for v in row['s']];u=read_function(row['stationary_scalar'])
        require(sum(s,ZERO)==ONE and all(not s[i] for i in range(3) if i not in I),'support/normalization identity')
        if sector=='neutral':diag=[v/2 for v in a];off=[(b[j]+c[j])/2 for j in range(3)];delta=[v/2 for v in d];weight=[-ONE]*3
        else:diag=[R[j][j] for j in range(3)];off=[R[i][j]+(1/alpha[k]+1/beta[k])/2 for k,(i,j) in enumerate(EDGES)];delta=[(1/alpha[j]-1/beta[j])/2 for j in range(3)];weight=w
        Q=[[ZERO for _ in range(3)] for _ in range(3)]
        for i in range(3):Q[i][i]=diag[i]
        for k,(i,j) in enumerate(EDGES):Q[i][j]=Q[j][i]=off[k]+(delta[k]*eps[i]*eps[j] if kind=='collinear' else ZERO)
        if kind=='triangle':
            correction=[delta[0]*delta[2]/delta[1],delta[0]*delta[1]/delta[2],delta[1]*delta[2]/delta[0]]
            for i in range(3):Q[i][i]-=correction[i]
        require(all(dot(Q[i],s)+weight[i]*u==ZERO for i in I),'stationary linear identities')
        reason=None;value=None
        if any(sign(s[i])<=0 for i in I):reason='nonpositive_support'
        if reason is None and kind=='collinear':cos=[rational(eps[i]*eps[j]) for i,j in EDGES]
        elif reason is None:
            g=[delta[0]*delta[2],delta[0]*delta[1],delta[1]*delta[2]];q=[sign(g[i])*g[i]*s[i] for i in range(3)]
            if any(sign(qi)<=0 for qi in q) or any(sign(sum(q,ZERO)-2*qi)<0 for qi in q):reason='unrealizable_triangle'
            else:
                sigma=sign(delta[0]*delta[1]*delta[2]);eta=[1,sigma*sign(delta[0]),sigma*sign(delta[2])]
                cos=[eta[i]*eta[j]*(q[3-i-j]**2-q[i]**2-q[j]**2)/(2*q[i]*q[j]) for i,j in EDGES]
        if reason is None and sector=='neutral':value=u
        elif reason is None:
            T=dot(w,s)+u*sum(w,ZERO)
            if not T:reason='zero_normalization'
            else:
                re=[dot(R[i],s) for i in range(3)];squared_norm=sum(((re[i]+u*w[i])**2 for i in range(3)),ZERO)
                for k,(i,j) in enumerate(EDGES):squared_norm+=s[i]*s[j]*((ONE+cos[k])/alpha[k]**2+(ONE-cos[k])/beta[k]**2)
                if sign(T*T-squared_norm)<0:reason='indefinite_gram'
                else:value=u/(2*T)
        require(row.get('accepted') is (reason is None) and row.get('rejection')==reason,'incorrect acceptance/rejection decision')
        if value is not None:
            require(read_function(row['value'])==value,'candidate value identity');v=limit(value)
            require(F(row['value_limit'])==v,'incorrect exact limit');values.append(v);counts[sector]=counts.get(sector,0)+1
        else:rejects[reason]=rejects.get(reason,0)+1
    require(bool(values),'no physical candidates')
    m=min(values);require(F(cert['minimum'])==m,'declared minimum is wrong')
    cls='positive' if m>0 else 'flat' if m==0 else 'negative';require(cert.get('classification')==cls,'wrong sign classification')
    return dict(verified=True,minimum=str(m),classification=cls,total_branches=28,accepted=counts,rejections=rejects,theorem_dependency='Full-domain completeness and exact-deformation lemmas; not an unconditional SOS certificate')

if __name__=='__main__':
    import json,sys
    if len(sys.argv)!=2:raise SystemExit('Usage: python verify_decision.py certificate.json')
    with open(sys.argv[1]) as stream:cert=json.load(stream)
    print(json.dumps(verify_decision(cert),indent=2))
