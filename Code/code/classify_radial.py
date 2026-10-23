"""Command-line entry point for the new full-domain original-norm formulation."""
import argparse,json,sys
from pathlib import Path
from full_radial_infinitesimal import run_radial

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path);parser.add_argument('--keep-functions',action='store_true');parser.add_argument('--regular-only',action='store_true',help='No deformation; refuse any singular branch system')
    args=parser.parse_args()
    try:
        p=json.loads(args.input.read_text());result=run_radial(**{k:p[k]for k in ('a','b','c','d')},deform=not args.regular_only,retain_functions=args.keep_functions);result['verification_boundary']='Research-theorem result, not an independent formal or human proof review.'
        text=json.dumps(result,default=str,indent=2)
        if args.output:args.output.write_text(text)
        print(json.dumps({k:result[k]for k in ('status','minimum','classification','systems','verification_boundary')},default=str,indent=2));return 0
    except (ValueError,TypeError,KeyError,ArithmeticError,OSError) as e:print(f'No classification: {e}',file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
