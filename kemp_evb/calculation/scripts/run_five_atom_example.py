#!/usr/bin/env python3
"""Run the reduced five-atom Kemp EVB teaching calculation."""
import argparse
from datetime import datetime
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from evb_teaching import run
from kemp_inputs import audit_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['Des27_7_IN'], default='Des27_7_IN', help='select the supplied Kemp reference case')
    parser.add_argument('--coupling', type=float, default=10.0, help='Illustrative H12, kcal/mol')
    parser.add_argument('--offset', type=float, default=0.0, help='Illustrative extra product offset, kcal/mol')
    parser.add_argument('--output', type=Path, help='New output directory; default includes date/time')
    parser.add_argument('--check', action='store_true', help='Check the bundled mapped inputs; do not calculate')
    args = parser.parse_args()
    if args.check:
        audit = audit_case(args.case)
        print(f"Input check passed: {audit['case']}; {audit['atoms']} atoms in source system; {audit['reactive_atoms']} mapped reactive atoms.")
        print('The live teaching calculation uses the five mapped atoms described in evb_teaching.py.')
        return
    stamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    output = args.output or HERE.parent / 'runs' / f'evb-demo-{stamp}'
    output = output if output.is_absolute() else (Path.cwd() / output).resolve()
    rows, summary = run(args.case, output, args.coupling, args.offset)
    print(f"Saved {len(rows)} points in {output}")
    print(f"Open {output / 'teaching_evb.html'}")
    print('This is a reduced teaching scan, not an enzyme free-energy profile.')


if __name__ == '__main__':
    main()
