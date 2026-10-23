"""Exact regular-first evaluation, with exact infinitesimal fallback.

No floating arithmetic or proximity threshold enters the classification. The
inherited solver is unchanged. Both routes rely on the manuscript's theorem;
this dispatcher is an efficiency wrapper, not independent proof verification.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from full_radial_infinitesimal import run_radial


def classify(parameters: dict, *, retain_functions: bool = False) -> dict:
    """Reject malformed inputs; fall back only on exact arithmetic singularities."""
    if set(parameters) != {'a', 'b', 'c', 'd'}: raise ValueError('Expected exactly the keys a, b, c, d.')
    t0 = time.perf_counter(); fallback_reason = None
    try:
        result = run_radial(**parameters, deform=False, retain_functions=retain_functions)
        route = 'exact_regular'
    except ArithmeticError as exc:
        fallback_reason = str(exc)
        result = run_radial(**parameters, deform=True, retain_functions=retain_functions)
        route = 'exact_infinitesimal_fallback'
    result.update(evaluation_route=route, fallback_reason=fallback_reason, elapsed_seconds=time.perf_counter()-t0,
                  arithmetic='exact rational / exact ordered rational-function arithmetic',
                  verification_boundary='Depends on the written universal theorem; not an independent proof audit.')
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path); parser.add_argument('--output', type=Path)
    parser.add_argument('--keep-functions', action='store_true'); args = parser.parse_args()
    try:
        result = classify(json.loads(args.input.read_text()), retain_functions=args.keep_functions)
        text = json.dumps(result, default=str, indent=2) + '\n'
        if args.output: args.output.write_text(text)
        print(json.dumps({k: result[k] for k in ('minimum','classification','evaluation_route','elapsed_seconds')},default=str,indent=2))
        return 0
    except (ValueError, TypeError, KeyError, ArithmeticError, OSError) as exc:
        print(f'No classification: {exc}', file=sys.stderr); return 2

if __name__ == '__main__': raise SystemExit(main())
