"""Classify an exact rational Weinberg quartic and emit a verified transcript.

Example: python classify_full.py input.json --certificate result.json
Input has four arrays a,b,c,d in edge order (12,23,31); use rational strings.
This implements the full-domain theorem draft, not a formal proof assistant.
"""
import argparse,json,sys
from pathlib import Path
from decision_certificate import export_decision
from verify_decision import verify_decision

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('--certificate',type=Path,required=True);args=parser.parse_args()
    try:
        p=json.loads(args.input.read_text());cert=export_decision(**p);result=verify_decision(cert);args.certificate.write_text(json.dumps(cert,indent=2));print(json.dumps(result,indent=2))
    except (ValueError,TypeError,ArithmeticError,KeyError,OSError) as exc:print(str(exc),file=sys.stderr);return 2
    return 0
if __name__=='__main__':raise SystemExit(main())
