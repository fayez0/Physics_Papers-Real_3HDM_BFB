"""Sharp original-field SOS identity for a nontrivial singular radial matrix."""
from pathlib import Path
from fractions import Fraction as F
import json
from field_sos import addto,bases,verify_sos,field_polynomial
ROOT=Path(__file__).parent;OUT=ROOT/'checkpoint04'

def build():
    p=dict(a=[2,5,10],b=[3,5,4],c=[-2,1,4],d=['1/3','2/5','3/7']);B=bases();lookup={mon:(k,i)for k,block in enumerate(B)for i,mon in enumerate(block)};grams=[[[F(0)for _ in block]for _ in block]for block in B];squares=[]
    def square(poly,w):
        active={lookup[m][0]for m in poly}
        if len(active)!=1:raise ValueError('Square mixes parity blocks')
        k=active.pop()
        for u,a in poly.items():
            for v,b in poly.items():grams[k][lookup[u][1]][lookup[v][1]]+=w*a*b
        squares.append(dict(weight=str(w),polynomial=[dict(monomial=list(m),coefficient=str(c))for m,c in poly.items()]))
    r=[{(4*i+j,4*i+j):F(1)for j in range(4)}for i in range(3)];rad={};addto(rad,r[0],F(4));addto(rad,r[1],F(-13));square(rad,F(1,102));square(r[2],F(212,51))
    for i,w in [(0,F(118,51)),(1,F(169,51))]:
        for a in range(4):
            for b in range(4):square({tuple(sorted((4*i+a,8+b))):F(1)},w)
    # det(Phi1,Phi2)=(x0+ix1)(x6+ix7)-(x2+ix3)(x4+ix5)
    dre={(0,6):F(1),(1,7):F(-1),(2,4):F(-1),(3,5):F(1)};dim={(0,7):F(1),(1,6):F(1),(2,5):F(-1),(3,4):F(-1)};square(dre,F(7,3));square(dim,F(7,3))
    for i,j,wre,wim in [(0,1,F(2,3),F(0)),(1,2,F(7,5),F(3,5)),(2,0,F(31,7),F(25,7))]:
        re={};im={}
        for k in (0,2):
            for u,v,sg in [(0,0,1),(1,1,1)]:addto(re,{tuple(sorted((4*i+k+u,4*j+k+v))):F(sg)})
            for u,v,sg in [(0,1,1),(1,0,-1)]:addto(im,{tuple(sorted((4*i+k+u,4*j+k+v))):F(sg)})
        if wre:square(re,wre)
        if wim:square(im,wim)
    cert=dict(schema='3hdm-original-fields-sos-v1',coefficients={k:list(map(str,v))for k,v in p.items()},lower_bound='43/51',sphere_multiplier_power=0,blocks=[dict(monomials=list(map(list,block)),gram=[[str(v)for v in row]for row in Q])for block,Q in zip(B,grams)],explicit_squares=squares)
    report=verify_sos(cert);(OUT/'sos_singular_A_sharp.json').write_text(json.dumps(cert,indent=2));(OUT/'singular_A_sharp_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':build()
