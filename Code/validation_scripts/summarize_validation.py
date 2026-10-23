"""Summarize saved rows; never infer an unexecuted test or an external result."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json, statistics
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'validation/validation_scripts'

def main():
    rows=[json.loads(s) for s in (D/'records.jsonl').read_text().splitlines()]
    expected=json.loads((D/'inputs.json').read_text())
    if len(rows)!=len(expected) or [x['id']for x in rows]!=[x['id']for x in expected]:raise ValueError('Incomplete or reordered final records.')
    for x in rows:
        if not x['completed']:raise ValueError(f'Uncompleted input: {x["id"]}')
        if x.get('predicted_shifted_minimum') and F(x['minimum'])!=F(x['predicted_shifted_minimum']):raise ValueError('Shift covariance failed.')
    cohorts={}
    for key in dict.fromkeys(x['cohort']for x in rows):
        sub=[x for x in rows if x['cohort']==key]
        cohorts[key]=dict(inputs=len(sub),max_abs_difference=max(abs(x['difference_absolute'])for x in sub),max_scaled_difference=max(abs(x['difference_scaled'])for x in sub),classification_counts=dict(Counter(x['classification']for x in sub)),branch_counts=dict(Counter(x['branch']for x in sub)))
    starts=[v for x in rows for v in x['starts']]
    routes={}
    for route in sorted(set(x['dispatch_route']for x in rows)):
        sub=[x for x in rows if x['dispatch_route']==route]
        routes[route]=dict(inputs=len(sub),median_dispatch_seconds=statistics.median(x['exact_dispatch_seconds']for x in sub),median_deformed_seconds=statistics.median(x['exact_deformed_seconds']for x in sub),median_paired_deformed_over_dispatch=statistics.median(x['exact_deformed_seconds']/x['exact_dispatch_seconds']for x in sub))
    summary=dict(inputs=len(rows),starts=len(starts),completed_inputs=len(rows),comparisons=dict(Counter(x['comparison']for x in rows)),cohorts=cohorts,
                 optimizer_non_success=sum(not x['success']for x in starts),optimizer_exceptions=sum(x.get('status')=='exception'for x in starts),
                 exact_route_agreements=len(rows),routes=routes,independent_sharp_reference_rows=sum(x['reference_kind']=='independent_sharp_analytic_minimum'for x in rows),
                 solver_centered_rows=sum(x['reference_kind']=='solver_centered_covariance_only'for x in rows),
                 classification_counts=dict(Counter(x['classification']for x in rows)),branch_counts=dict(Counter(x['branch']for x in rows)),
                 numeric_is_global_certificate=False,reviewer_report_rows_used=0,timing_caveat='In-process, imports excluded, four concurrent workers; not a matched StableWein runtime comparison.')
    (D/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (D/'stablewein_full_inputs.json').write_text(json.dumps([{k:x[k]for k in ('id','couplings','minimum','cohort','reference_kind')}for x in rows],indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
