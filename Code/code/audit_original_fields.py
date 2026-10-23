"""Independent symbolic coefficient audit and exact SOS/separator verification."""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import json,time
import sympy as sp
from field_sos import field_polynomial,verify_sos,verify_separator,exact_field_value
ROOT=Path(__file__).parent;OUT=ROOT/'checkpoint04'

def run():
    tic=time.monotonic();x=sp.symbols('x0:12',real=True);Phi=[sp.Matrix([x[4*i]+sp.I*x[4*i+1],x[4*i+2]+sp.I*x[4*i+3]])for i in range(3)];K=sp.Matrix(3,3,lambda i,j:(Phi[i].conjugate().T*Phi[j])[0]);rad=[sp.expand(K[i,i])for i in range(3)];edges=[(0,1),(1,2),(2,0)];symterms=[r*r/2 for r in rad]+[rad[i]*rad[j]for i,j in edges]+[sp.expand(K[i,j]*sp.conjugate(K[i,j])) for i,j in edges]+[sp.expand((K[i,j]**2+sp.conjugate(K[i,j])**2)/2)for i,j in edges]+[-sum(rad)**2];identities=[]
    for k,target in enumerate(symterms):
        coeff=[[0]*3 for _ in range(4)]
        if k<12:coeff[k//3][k%3]=1
        output=field_polynomial(*coeff,lower=int(k==12));poly=sp.Poly(sp.expand(target),*x);expected={tuple(i for i,e in enumerate(exponent)for _ in range(e)):F(c) for exponent,c in poly.terms() if c}
        if output!=expected:raise ArithmeticError('Original potential coefficient mismatch: '+str(k))
        identities.append(dict(basis_index=k,passed=True,nonzero_monomials=len(expected)))
    verified=[]
    for path in sorted((OUT/'sos_certificates').glob('*.json'))+sorted(OUT.glob('sos_symmetric*.json'))+sorted(OUT.glob('sos_singular*.json')):
        c=json.loads(path.read_text())
        if c.get('schema')!='3hdm-original-fields-sos-v1':continue
        result=verify_sos(c);verified.append(dict(file=str(path.relative_to(ROOT)),**result))
    sep=json.loads((OUT/'quartic_sos_obstruction.json').read_text());separator=verify_separator(sep);mutations=[]
    base=json.loads(next((OUT/'sos_certificates').glob('*.json')).read_text())
    def trial(name,mutate,source=base,fn=verify_sos):
        c=deepcopy(source);mutate(c)
        try:fn(c)
        except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError):mutations.append(dict(name=name,rejected=True));return
        raise ArithmeticError('Corrupted certificate accepted: '+name)
    trial('changed_lower_bound',lambda c:c.__setitem__('lower_bound',str(F(c['lower_bound'])+1)))
    trial('omitted_block',lambda c:c['blocks'].pop())
    trial('wrong_monomials',lambda c:c['blocks'][0]['monomials'].__setitem__(0,[0,0]))
    trial('negative_diagonal',lambda c:c['blocks'][0]['gram'][0].__setitem__(0,'-1'))
    trial('nonsymmetric_gram',lambda c:c['blocks'][0]['gram'][0].__setitem__(1,'999'))
    trial('changed_potential',lambda c:c['coefficients']['a'].__setitem__(0,str(F(c['coefficients']['a'][0])+1)))
    trial('invalid_multiplier',lambda c:c.__setitem__('sphere_multiplier_power',-1))
    trial('separator_missing_moment',lambda c:c['moment_values'].pop(),source=sep,fn=verify_separator)
    trial('separator_false_pairing',lambda c:c.__setitem__('separation','-1'),source=sep,fn=verify_separator)
    output=dict(status='passed',symbolic_original_potential_identities=identities,sos_certificates=verified,separator=separator,mutation_tests=mutations,assertions_enabled=__debug__,elapsed_seconds=time.monotonic()-tic);(OUT/('original_fields_optimized.json' if not __debug__ else 'original_fields_normal.json')).write_text(json.dumps(output,indent=2));print(json.dumps({k:v for k,v in output.items() if k not in ('sos_certificates','symbolic_original_potential_identities')},indent=2))
if __name__=='__main__':run()
