"""Direct symbolic attainment checks in actual complex doublets."""
from pathlib import Path
import json,sympy as s
ROOT=Path(__file__).parent;OUT=ROOT/'checkpoint04'

def value(X,a,b,c,d):
    K=s.simplify(X.conjugate().T*X);r=[K[i,i]for i in range(3)];V=sum(a[i]*r[i]**2/2 for i in range(3))
    for e,(i,j) in enumerate([(0,1),(1,2),(2,0)]):V+=b[e]*r[i]*r[j]+c[e]*K[i,j]*s.conjugate(K[i,j])+d[e]*(K[i,j]**2+s.conjugate(K[i,j])**2)/2
    return s.simplify(V/s.trace(K)**2),K

def run():
    X=s.Matrix([[0,0,0],[s.sqrt(13),2*s.I,0]]);val,K=value(X,[2,5,10],[3,5,4],[-2,1,4],[s.Rational(1,3),s.Rational(2,5),s.Rational(3,7)])
    if val!=s.Rational(43,51) or s.det(s.Matrix([[2,3,4],[3,5,5],[4,5,10]]))!=0:raise ArithmeticError('Singular-A attainment failed')
    w=-s.Rational(1,2)+s.I*s.sqrt(3)/2;X=s.Matrix([[1/s.sqrt(5),w/s.sqrt(5),w**2/s.sqrt(5)],[s.sqrt(s.Rational(2,15))]*3]);val2,K2=value(X,[3]*3,[-2]*3,[1]*3,[s.Rational(1,2)]*3)
    if val2!=-s.Rational(7,60) or s.simplify(s.trace(K2*K2))!=s.Rational(13,25) or s.simplify(s.det(K2))!=0:raise ArithmeticError('Complex charged attainment failed')
    result=dict(status='passed',singular_A_minimum=str(val),singular_A_witness_trace=str(s.trace(K)),complex_charged_minimum=str(val2),complex_charged_Gram_trace=str(s.trace(K2)),complex_charged_Gram_square_trace='13/25');(OUT/'sharp_field_attainment.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':run()
