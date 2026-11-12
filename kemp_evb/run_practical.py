#!/usr/bin/env python3
"""Open one of the saved Kemp EVB workshop examples."""
import argparse
import hashlib
import json
from pathlib import Path
import webbrowser

ROOT = Path(__file__).resolve().parent
PAGES = {
    "core": "viewers/evb_core.html",
    "profile": "viewers/source_profile.html",
    "scan": "viewers/parameter_scan.html",
}


def verify_files():
    manifest = json.loads((ROOT / "workshop_manifest.json").read_text(encoding="utf-8"))
    failures = []
    for name, digest in manifest["files"].items():
        path = ROOT / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            failures.append(name)
    if failures:
        raise SystemExit("Workshop files missing or changed: " + ", ".join(failures))
    print(f"Kemp workshop files verified: {len(manifest['files'])}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", choices=(*PAGES, "check"), default="check")
    args = parser.parse_args()
    verify_files()
    if args.step != "check":
        page = ROOT / PAGES[args.step]
        print(f"Opening saved Kemp example: {page}")
        webbrowser.open(page.as_uri())
    else:
        print("Open index.html for the hands-on. This launcher checks files; it does not run simulation stages.")


if __name__ == "__main__":
    main()
