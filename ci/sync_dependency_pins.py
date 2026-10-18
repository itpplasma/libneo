#!/usr/bin/env python3
"""Regenerate Nix dependency pins from fpm.toml, or check their consistency."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import tomllib

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
dependencies = tomllib.loads((root / "fpm.toml").read_text())["dependencies"]
flake = root / "flake.nix"
text = flake.read_text()
updated = text
for name in ("fortio", "fortnum"):
    revision = dependencies[name]["rev"]
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise SystemExit(f"{name}: expected an exact commit in fpm.toml")
    updated, count = re.subn(
        rf'github:lazy-fortran/{name}/[^"\s]+',
        f"github:lazy-fortran/{name}/{revision}", updated,
    )
    if count != 1:
        raise SystemExit(f"{name}: expected one Nix input URL")

if args.check:
    if text != updated:
        raise SystemExit("Nix inputs differ from fpm.toml; run ci/sync_dependency_pins.py")
else:
    flake.write_text(updated)
    subprocess.run([
        "nix", "flake", "update", "fortio", "fortnum",
    ], cwd=root, check=True)

lock = json.loads((root / "flake.lock").read_text())
for name in ("fortio", "fortnum"):
    node = lock["nodes"][name]
    for kind in ("original", "locked"):
        if node[kind]["rev"] != dependencies[name]["rev"]:
            raise SystemExit(f"{name}: Nix {kind} revision differs from fpm.toml")
print("Fortio and Fortnum pins agree across the build entrypoints.")
