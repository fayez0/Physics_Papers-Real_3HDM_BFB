"""Independent lower/upper certificates, NOT a proof of universal completeness."""
from pathlib import Path
from fractions import Fraction as F
from field_sos import GramSearch,verify_sos,find_upper_witness,exact_field_value,verify_separator
from exact_regular import enumerate_regular
from full_infinitesimal import full_minimum
import json,random,time,sys
ROOT=Path(__file__).parent; OUT=ROOT/'checkpoint04'; OUT.mkdir(exist_ok=True)

def dataset():
    external=json.loads((ROOT/'checkpoint03/external_benchmark_inputs.json').read_text()); rows=[dict(id=row['id'],couplings=row['couplings'],source='checkpoint03_selected_charged_pair') for row in external]
    rng=random.Random(2026101612)
    for i in range(16): rows.append(dict(id=f'random_{i:02d}',couplings=dict(a=[str(F(rng.randint(1,31),5)) for _ in range(3)],b=[str(F(rng.randint(-28,24),7)) for _ in range(3)],c=[str(F(rng.randint(-25,29),9)) for _ in range(3)],d=[str(F(rng.choice([-1,1])*rng.randint(1,23),11)) for _ in range(3)]),source='seeded_unshifted_rational_box'))
    return rows

def run(start=0,stop=28):
    certdir=OUT/'sos_certificates';certdir.mkdir(exist_ok=True); results=[]; solvers=[GramSearch(),GramSearch(1)];starttime=time.monotonic()
    for index,row in list(enumerate(dataset()))[start:stop]:
        t=time.monotonic();p=row['couplings']; reference=enumerate_regular(**p)
        if reference['minimum'] is None:
            exact=full_minimum(**p);m=exact['minimum'];branch=exact['best_branch_limit']; route='infinitesimal'
        else: m=reference['minimum']['value'];branch=reference['minimum'];route='regular'
        margin=min(F(1,1000),abs(m)/4) if m else F(1,1000); lower=m-margin; attempts=[];certificate=None
        for power,solver in enumerate(solvers):
            result=solver.search(p,lower,maxiter=3000 if power==0 else 6000);attempts.append(dict(multiplier_power=power,status=result['status'],iterations=result['iterations'],final_search_diagnostic=result['log'][-1]))
            if result['status']=='exact_sos_certificate':certificate=result['certificate'];break
        upper=find_upper_witness(p,starts=6,seed=2026101610+index); U=F(upper['upper_bound']); entry=dict(index=index,**row,predicted_minimum=str(m),prediction_route=route,predicted_sector=branch['sector'],predicted_branch=branch['branch'],attempts=attempts,upper_witness=upper)
        if U<m: raise ArithmeticError(f'Independent field below prediction: {index}')
        if certificate:
            check=verify_sos(certificate); name=f'{index:03d}.json'; (certdir/name).write_text(json.dumps(certificate,indent=2));entry.update(lower_certificate=name,certificate_check=check,independent_interval=[str(lower),str(U)],interval_width=str(U-lower),sign_certified=(lower>0 or U<0))
        else:entry.update(sign_certified=U<0,lower_bound_status='no_certificate_found_no_nonnegativity_conclusion')
        entry['elapsed_seconds']=time.monotonic()-t;results.append(entry);print(index,branch['sector'],branch['branch'],str(m),'cert',bool(certificate),'power',power,'secs',round(entry['elapsed_seconds'],3),flush=True)
        (OUT/f'independent_sos_{start}_{stop}.json').write_text(json.dumps(dict(status='completed_current_rows',completed=len(results),requested=stop-start,rows=results,elapsed_seconds=time.monotonic()-starttime),indent=2))
    return results
if __name__=='__main__': run(*list(map(int,sys.argv[1:3])))
