#!/usr/bin/env python3
"""Open the offline workshop, or verify its delivered files (standard library only)."""
import argparse
import hashlib
import json
from pathlib import Path
import webbrowser

ROOT=Path(__file__).resolve().parent


def check_assets():
    manifest=json.loads((ROOT/'workshop_manifest.json').read_text())
    failures=[]
    for name,digest in manifest['files'].items():
        path=ROOT/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:failures.append(name)
    if failures:raise ValueError('Missing/changed workshop files: '+', '.join(failures)+'. Pull the complete updated folder.')
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true')
    args=parser.parse_args();manifest=check_assets()
    print(f'Workshop files verified: {len(manifest["files"])}. Kemp EVB hands-on.')
    if not args.check:
        page=ROOT/'index.html';print(f'Open: {page}')
        webbrowser.open(page.as_uri())

if __name__=='__main__':
    try:main()
    except (ValueError,FileNotFoundError) as error:raise SystemExit(str(error))
