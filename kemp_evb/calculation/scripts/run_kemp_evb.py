#!/usr/bin/env python3
"""Full-enzyme Kemp EVBOpenMM adaptation: relaxation through mapping and profile.
Default: complete short full-enzyme run. This is not exact Q6 reproduction.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]/'full_enzyme'
sys.path.insert(0,str(ROOT/'scripts'))

def main():
    if '--check' in sys.argv:
        sys.argv.remove('--check');sys.argv.extend(['--mode','validate'])
    def supplied(flag):return any(x==flag or x.startswith(flag+'=') for x in sys.argv[1:])
    if not supplied('--mode') and not any(x in sys.argv for x in ('-h','--help')):
        sys.argv.extend(['--mode','run'])
    if not supplied('--config'):
        sys.argv.extend(['--config',str(ROOT/'configs/full_enzyme_short.json')])
    from full_enzyme import main as run
    return run()

if __name__=='__main__':
    try:main()
    except (ValueError,FileNotFoundError,FileExistsError) as e:raise SystemExit(str(e))
