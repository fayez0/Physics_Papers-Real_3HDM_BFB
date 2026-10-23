"""Exact attaining Gram matrices, including singular limits, for the real Weinberg quartic.

The generator uses the checkpoint-04 full theorem to select a branch. Its output
is an exact feasible-field witness: it proves attainment of the reported value,
not a lower bound independently of the full-domain theorem. The separate
verify_minimizer.py does not import a branch solver or use floating arithmetic.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
import sympy as sp
from full_radial_infinitesimal import run_radial,RF,E,Z,O,EDGES,P0,q,sign,limit,angular_inverse

def rf_expr(text: str):
    """Parse only internal expressions produced by the retained branch solver."""
    return RF.from_expr(sp.sympify(text))

def signed_sqrt_sum(a: int, x, b: int, y) -> int:
    """Eventual sign of a sqrt(x)+b sqrt(y), with nonnegative x,y and a,b in {-1,0,1}."""
    sx,sy=sign(x),sign(y)
    if sx<0 or sy<0: raise ArithmeticError('Negative squared half-phase amplitude')
    a=a if sx else 0; b=b if sy else 0
    if not a:return b
    if not b:return a
    return a if a==b else a*sign(x-y)

def phase_data(weights, edge_phase, kind: str, epsilon=None):
    """Canonical doubled phases z_i; return Re(z_i), signs of sin(theta_i), and cos(2Delta)."""
    if kind=='collinear':
        re=[RF(x) for x in epsilon];sigma=[0 if x==1 else 1 for x in epsilon]
    else:
        g=[abs_rf(edge_phase[0]*edge_phase[2]),abs_rf(edge_phase[0]*edge_phase[1]),abs_rf(edge_phase[1]*edge_phase[2])]
        Q=[g[i]*weights[i] for i in range(3)];Q1,Q2,Q3=Q
        if any(sign(x)<=0 for x in Q):raise ArithmeticError('Triangle needs full positive support')
        h=4*Q1*Q1*Q2*Q2-(Q3*Q3-Q1*Q1-Q2*Q2)**2
        if sign(h)<0:raise ArithmeticError('Impossible phase triangle')
        sg=sign(edge_phase[0]*edge_phase[1]*edge_phase[2]);eta=[1,sg*sign(edge_phase[0]),sg*sign(edge_phase[2])]
        re0=[O,(Q3**2-Q1**2-Q2**2)/(2*Q1*Q2),(Q2**2-Q1**2-Q3**2)/(2*Q1*Q3)]
        im0=[Z,1/(2*Q1*Q2),-1/(2*Q1*Q3)]
        re=[eta[i]*re0[i] for i in range(3)]
        sigma=[(sign(eta[i]*im0[i]) if h else (1 if re[i]==-O else 0)) for i in range(3)]
        # Nonzero h implies sqrt(h)>0 on a small positive interval; at z=-1 choose theta=pi/2.
        for i in range(3):
            if sign(O-re[i])<0 or sign(O+re[i])<0:raise ArithmeticError('Invalid doubled-phase real part')
        cos=[re[i]*re[j]+eta[i]*eta[j]*im0[i]*im0[j]*h for i,j in EDGES]
        return re,sigma,cos
    return re,sigma,[RF(epsilon[i]*epsilon[j]) for i,j in EDGES]

def abs_rf(x):return x if sign(x)>=0 else -x

def half_angle_signs(re,sigma,i: int,j: int):
    """Signs of cos(theta_i-theta_j), sin(theta_i-theta_j), without taking radicals."""
    cr=signed_sqrt_sum(1,(O+re[i])*(O+re[j]),sigma[i]*sigma[j],(O-re[i])*(O-re[j]))
    si=signed_sqrt_sum(sigma[i],(O-re[i])*(O+re[j]),-sigma[j],(O+re[i])*(O-re[j]))
    return cr,si

def root_record(square,sgn: int):
    if sign(square)<0:raise ArithmeticError('Negative square before limit')
    value=limit(square)
    if value<0:raise ArithmeticError('Negative square limit')
    if value and sgn==0:raise ArithmeticError('Nonzero amplitude has zero sign')
    return {'sign':sgn if value else 0,'square':str(value)}

def generate(p: dict, branch_index=None, *, _result=None):
    if set(p)!=set(('a','b','c','d')):raise ValueError('Input keys must be exactly a,b,c,d')
    if any(len(v)!=3 for v in p.values()):raise ValueError('Each coefficient array must have length three')
    p={k:[str(F(v)) if not isinstance(v,(float,bool)) else reject_float() for v in p[k]] for k in ('a','b','c','d')}
    result=run_radial(**p,retain_functions=True) if _result is None else _result
    row=result['best'] if branch_index is None else result['records'][branch_index]
    if not row.get('accepted'):raise ValueError('Requested branch is not eventually physical')
    branch_number=next(i for i,b in enumerate(result['records']) if b is row)
    pp={k:[q(v)+E*q(P0[k][i]) for i,v in enumerate(p[k])] for k in p};c,d=pp['c'],pp['d']
    r=[rf_expr(x) for x in row['radii']];target=limit(rf_expr(row['value_function']))
    zero={'real':{'sign':0,'square':'0'},'imag':{'sign':0,'square':'0'}}
    entries=[json.loads(json.dumps(zero)) for _ in EDGES]
    label=row['label'];kernel_s=None
    if label=='charged_orthogonal_pair':pass
    elif label=='charged_split':
        pair=tuple(row['pair']);e=EDGES.index(pair);entries[e]['real' if row['phase_sign']==1 else 'imag']=root_record(r[pair[0]]*r[pair[1]],1)
    else:
        if label.startswith('neutral'):
            s=r;nu=O;alpha=[O]*3;beta=[O]*3;delta=d
            kind='triangle' if label.endswith('triangle') else 'collinear';eps=row.get('phase_signs')
        else:
            kind='triangle' if label.endswith('triangle') else 'collinear';eps=row.get('phase_signs');G=angular_inverse(c,d,kind,eps)
            z=[sum((G[i][j]*r[j] for j in range(3)),Z) for i in range(3)];total=sum(z,Z)
            if not total:raise ArithmeticError('Zero kernel normalization')
            s=[v/total for v in z];kernel_s=[str(limit(x)) for x in s];nu=-total
            alpha=[c[e]+d[e] for e in range(3)];beta=[c[e]-d[e] for e in range(3)];delta=[(1/alpha[e]-1/beta[e])/2 for e in range(3)]
        re,sigma,cos=phase_data(s,delta,kind,eps)
        for e,(i,j) in enumerate(EDGES):
            cr,si=half_angle_signs(re,sigma,i,j);rs=nu**2*s[i]*s[j]*(O+cos[e])/(2*alpha[e]**2);ims=nu**2*s[i]*s[j]*(O-cos[e])/(2*beta[e]**2)
            entries[e]={'real':root_record(rs,cr*sign(nu/alpha[e])),'imag':root_record(ims,si*sign(nu/beta[e]))}
    return {'schema':'weinberg-attaining-gram-v1','coefficients':p,'edge_order':[[1,2],[2,3],[3,1]],'claimed_value':str(target),'reported_global_minimum':str(result['minimum']),'is_minimum_value_branch':target==result['minimum'],'source_branch_index':branch_number,'source_branch':label,'diagonal':[str(limit(x)) for x in r],'off_diagonal':entries,'limiting_kernel_norms':kernel_s,'verification_boundary':'An exact feasible Gram matrix proves attainment, not a universal lower bound. Global optimality uses the full-domain theorem unless a separate lower certificate is supplied.'}

def reject_float():raise TypeError('Use integers or exact strings, not floats or booleans')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--branch',type=int);args=parser.parse_args()
    try:
        cert=generate(json.loads(args.input.read_text()),args.branch);args.output.write_text(json.dumps(cert,indent=2)+'\n');print(json.dumps({k:cert[k] for k in ('claimed_value','source_branch','is_minimum_value_branch')},indent=2));return 0
    except (ValueError,TypeError,KeyError,ArithmeticError,OSError) as exc:print('No witness:',exc);return 2
if __name__=='__main__':raise SystemExit(main())
