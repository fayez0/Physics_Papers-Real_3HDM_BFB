"""Check proposed value-function degree bounds along exact linear deformations."""
import json,time
from pathlib import Path
from fractions import Fraction as F
import sympy as sp
from full_infinitesimal import full_minimum, RF
ROOT=Path(__file__).parent

def run():
    t=time.monotonic(); rows=[]; maxima={'neutral':[0,0],'charged':[0,0]}
    special=json.loads((ROOT/'checkpoint02'/'special_limits.json').read_text()); cases=[r['coefficients'] for r in special['rows'][::5]]+[dict(a=[3]*3,b=[-2]*3,c=[1]*3,d=['1/2']*3),dict(a=[0]*3,b=[0]*3,c=[0]*3,d=[0]*3)]
    for i,p in enumerate(cases):
        result=full_minimum(**p,keep_functions=True)
        for c in result['candidates']:
            f=RF.from_expr(sp.sympify(c['value_function'])); n=max([key[0] for key in f.numer.to_dict()] or [0]); d=max([key[0] for key in f.denom.to_dict()] or [0]); limits=(12,11) if c['sector']=='neutral' else (36,35)
            assert n<=limits[0] and d<=limits[1]; maxima[c['sector']]=[max(n,maxima[c['sector']][0]),max(d,maxima[c['sector']][1])]; rows.append(dict(case=i,sector=c['sector'],branch=c['branch'],numerator_degree=n,denominator_degree=d))
    report=dict(status='passed',input_count=len(cases),candidate_functions=len(rows),proved_upper_bounds={'neutral':[12,11],'charged':[36,35]},observed_maxima=maxima,rows=rows,elapsed_seconds=time.monotonic()-t); (ROOT/'checkpoint02'/'degree_bounds.json').write_text(json.dumps(report,indent=2)); print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
if __name__=='__main__': run()
