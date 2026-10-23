"""Exact original-field attainment verifier, independent of BFB solvers and SymPy.

Input entries are signed square roots of nonnegative rationals, never executable
code. Radical arithmetic uses rational-square tests, not approximate equality.
This proves a feasible field value and stationarity, NOT a global lower bound.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
from radical_arithmetic import Radical as R
EDGES=((0,1),(1,2),(2,0))

def require(test,message):
    if not test:raise ValueError(message)

def rational(x):
    if isinstance(x,(float,bool)):raise ValueError('Inexact or boolean rational')
    return F(x)

def root(record):
    require(set(record)=={'sign','square'},'Malformed radical');s=record['sign'];v=rational(record['square'])
    require(type(s) is int and s in (-1,0,1),'Invalid radical sign');require(v>=0,'Negative radical');require((v==0)==(s==0),'Noncanonical zero')
    return s*R.sqrt(v)

def gram(X):return [[sum((row[i].conjugate()*row[j] for row in X),R()) for j in range(3)] for i in range(3)]

def verify(cert, export_fields=True):
    require(cert.get('schema')=='weinberg-attaining-gram-v1','Unknown schema');require(cert.get('edge_order')==[[1,2],[2,3],[3,1]],'Wrong edge order')
    p=cert['coefficients'];require(set(p)==set(('a','b','c','d')),'Wrong coefficient set');require(all(len(v)==3 for v in p.values()),'Wrong coefficient length')
    a,b,c,d=([rational(v) for v in p[k]] for k in ('a','b','c','d'));r=[rational(v) for v in cert['diagonal']]
    require(len(r)==3 and len(cert['off_diagonal'])==3,'Wrong matrix dimensions');require(sum(r)==1,'Gram trace differs from one');require(all(v>=0 for v in r),'Negative doublet norm')
    K=[[R(r[i]) if i==j else R() for j in range(3)] for i in range(3)]
    for e,(i,j) in enumerate(EDGES):
        item=cert['off_diagonal'][e];require(set(item)=={'real','imag'},'Malformed matrix entry');K[i][j]=root(item['real'])+R.imaginary()*root(item['imag']);K[j][i]=K[i][j].conjugate()
    pairs=[(r[i]*r[j]-K[i][j]*K[j][i]).rational() for i,j in EDGES];require(all(x>=0 for x in pairs),'Negative principal pair minor')
    determinant=K[0][0]*(K[1][1]*K[2][2]-K[1][2]*K[2][1])-K[0][1]*(K[1][0]*K[2][2]-K[1][2]*K[2][0])+K[0][2]*(K[1][0]*K[2][1]-K[1][1]*K[2][0]);require(not determinant,'Gram determinant is nonzero')
    actual=R(sum(a[i]*r[i]**2/2 for i in range(3)))
    for e,(i,j) in enumerate(EDGES):actual+=b[e]*r[i]*r[j]+c[e]*K[i][j]*K[j][i]+d[e]*(K[i][j]*K[i][j]+K[j][i]*K[j][i])/2
    actual=actual.rational();claimed=rational(cert['claimed_value']);require(actual==claimed,'Original potential differs from claimed value')
    i=next(i for i,x in enumerate(r) if x>0);v=[K[i][j]*R.sqrt(1/r[i]) for j in range(3)];G=gram([v]);rem=[[K[i][j]-G[i][j] for j in range(3)] for i in range(3)];rd=[rem[i][i].rational() for i in range(3)];positive=[j for j,x in enumerate(rd) if x>0]
    if positive:
        k=positive[0];w=[rem[k][j]*R.sqrt(1/rd[k]) for j in range(3)];X=[v,w];rank=2
    else:
        require(all(not x for row in rem for x in row),'Unresolved residual');X=[v,[R(),R(),R()]];rank=1
    G=gram(X);require(all(K[i][j]==G[i][j] for i in range(3) for j in range(3)),'Factorized fields do not reproduce Gram')
    L=[[R() for _ in range(3)] for _ in range(3)]
    for i in range(3):L[i][i]=R(a[i]*r[i]+sum(b[e]*r[j] for e,(u,vv) in enumerate(EDGES) for j in [vv if u==i else u] if u==i or vv==i))
    for e,(i,j) in enumerate(EDGES):L[i][j]=c[e]*K[i][j]+d[e]*K[j][i];L[j][i]=L[i][j].conjugate()
    stationary=all(not (sum((X[u][k]*L[k][j] for k in range(3)),R())-2*actual*X[u][j]) for u in range(2) for j in range(3))
    out={'status':'exact_attainment_verified','value':str(actual),'rank':rank,'doublet_norms':[str(x) for x in r],'principal_pair_minors':[str(x) for x in pairs],'field_stationarity':stationary,'factorization_verified':True,'independent_global_lower_bound':False}
    if export_fields:out.update(fields=[[x.text() for x in row] for row in X],field_exact_terms=[[x.payload() for x in row] for row in X])
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('certificate',type=Path);parser.add_argument('--output',type=Path);args=parser.parse_args()
    try:
        result=verify(json.loads(args.certificate.read_text()));text=json.dumps(result,indent=2)
        if args.output:args.output.write_text(text+'\n')
        print(text);return 0
    except (ValueError,TypeError,KeyError,OSError,ArithmeticError) as exc:print('Certificate rejected:',exc);return 2
if __name__=='__main__':raise SystemExit(main())
