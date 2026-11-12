#!/usr/bin/env python3
"""Recompute the short Kemp exercises using bundled energies, not new enzyme MD."""
import argparse
import hashlib
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import webbrowser

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
SCRIPTS=REPO/'Kemp_EVB/scripts'


def check_compute():
    from start_workshop import check_assets
    check_assets()
    if not SCRIPTS.is_dir():raise ValueError('Recomputation needs the full EVB_MHETase checkout. The offline index/viewers work with this folder alone.')
    missing=[name for name in ('openmm','numpy','evbopenmm') if importlib.util.find_spec(name) is None]
    if missing:raise ValueError('Activate kemp-openmm and complete SETUP.md. Missing: '+', '.join(missing)+'. Saved viewers remain usable without these packages.')
    result=subprocess.run([sys.executable,'-c','import openmm,evbopenmm; print("OpenMM",openmm.version.version); print("EVBOpenMM",evbopenmm.__file__)'],check=True,capture_output=True,text=True)
    print(result.stdout.strip())


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step',choices=('core','profile','scan'),default='scan')
    parser.add_argument('--check',action='store_true');parser.add_argument('--no-open',action='store_true')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args(argv);check_compute()
    if args.check:return 0
    if args.step=='core':
        if args.output_dir:out=args.output_dir.resolve()
        else:
            index=1
            while (ROOT/'outputs'/f'core_{index:03d}').exists():index+=1
            out=ROOT/'outputs'/f'core_{index:03d}'
        if out.exists():raise ValueError('Core output already exists; choose a new --output-dir. Previous results were not changed.')
        command=[sys.executable,str(SCRIPTS/'kemp_protocol.py'),'--mode','evb-demo','--output-dir',str(out)]
        subprocess.run(command,cwd=REPO/'Kemp_EVB',check=True)
        if not args.no_open:webbrowser.open((out/'teaching_evb.html').as_uri())
    else:
        out=args.output_dir.resolve() if args.output_dir else ROOT/'outputs/analysis'
        out.mkdir(parents=True,exist_ok=True)
        source=ROOT/'data/short_run/replica_101_raw.json';target=out/source.name
        digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
        if target.exists() and digest(target)!=digest(source):raise ValueError('Output contains a different energy dataset; choose a new --output-dir. No records were overwritten.')
        if not target.exists():shutil.copy2(source,target)
        command=[sys.executable,str(SCRIPTS/'full_enzyme.py'),'--mode',args.step,'--output-dir',str(out),'--bins','18']
        if not args.no_open:command.append('--open')
        subprocess.run(command,cwd=REPO/'Kemp_EVB',check=True)
    print(f'Exercise outputs: {out}')
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,FileNotFoundError,subprocess.CalledProcessError) as error:raise SystemExit(f'Exercise stopped: {error}')
