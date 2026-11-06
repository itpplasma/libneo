#!/usr/bin/env python3
"""Regenerate Nix dependency pins from fpm.toml, or check their consistency.

fpm.toml pins Fortio, Fortnum and FortFEM. FortSparse and FortPlot follow the
revisions that FortFEM's CMake build fetches at its pinned commit.
"""

import argparse
import json
from pathlib import Path
import re
import subprocess
import tomllib

DIRECT = ("fortio", "fortnum", "fortfem")
TRANSITIVE = ("fortsparse", "fortplot")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
dependencies = tomllib.loads((root / "fpm.toml").read_text())["dependencies"]
flake = root / "flake.nix"
text = flake.read_text()


def nix_input_revision(name, source):
    match = re.findall(rf'github:lazy-fortran/{name}/([0-9a-f]{{40}})', source)
    if len(match) != 1:
        raise SystemExit(f"{name}: expected one Nix input URL")
    return match[0]


def fortfem_cmake_pins(revision):
    prefetch = subprocess.run(
        ["nix", "flake", "prefetch", "--json",
         f"github:lazy-fortran/fortfem/{revision}"],
        check=True, capture_output=True, text=True,
    )
    cmake = Path(json.loads(prefetch.stdout)["storePath"]) / "CMakeLists.txt"
    source = cmake.read_text()
    pins = {}
    for name in TRANSITIVE:
        match = re.search(
            rf"FetchContent_Declare\({name}\s+GIT_REPOSITORY\s+\S+\s+"
            rf"GIT_TAG\s+([0-9a-f]{{40}})", source)
        if match is None:
            raise SystemExit(f"{name}: no pinned GIT_TAG in FortFEM's CMakeLists.txt")
        pins[name] = match.group(1)
    return pins


expected = {}
for name in DIRECT:
    revision = dependencies[name]["rev"]
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise SystemExit(f"{name}: expected an exact commit in fpm.toml")
    expected[name] = revision
if args.check:
    for name in TRANSITIVE:
        expected[name] = nix_input_revision(name, text)
else:
    expected.update(fortfem_cmake_pins(expected["fortfem"]))

updated = text
for name, revision in expected.items():
    nix_input_revision(name, updated)
    updated = re.sub(
        rf'github:lazy-fortran/{name}/[^"\s]+',
        f"github:lazy-fortran/{name}/{revision}", updated,
    )

if args.check:
    if text != updated:
        raise SystemExit("Nix inputs differ from fpm.toml; run ci/sync_dependency_pins.py")
else:
    flake.write_text(updated)
    subprocess.run(["nix", "flake", "update", *expected], cwd=root, check=True)

lock = json.loads((root / "flake.lock").read_text())
for name, revision in expected.items():
    node = lock["nodes"][name]
    for kind in ("original", "locked"):
        if node[kind]["rev"] != revision:
            raise SystemExit(f"{name}: Nix {kind} revision differs from the pins")
print("Fortio, Fortnum and FortFEM pins agree across the build entrypoints.")
