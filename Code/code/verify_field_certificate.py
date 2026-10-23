"""Check an original-field SOS or separating-moment certificate exactly."""
import argparse,json,sys
from pathlib import Path
from field_sos import verify_sos,verify_separator

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('certificate',type=Path);args=parser.parse_args()
    try:
        cert=json.loads(args.certificate.read_text());result=(verify_separator if cert.get('schema')=='3hdm-sos-separator-v1' else verify_sos)(cert);print(json.dumps(result,indent=2));return 0
    except (ValueError,TypeError,KeyError,ArithmeticError,OSError,IndexError)as e:print(f'Invalid certificate: {e}',file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
