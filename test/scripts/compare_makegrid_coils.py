#!/usr/bin/env python3
"""Cross-check the Fortran MAKEGRID reader against python/libneo/coils.py.

A synthetic coils file with groups, comments, two filament blocks, an
unterminated tail coil and a repeated closing point is parsed by both readers.
Coil geometry and currents must agree bitwise, group ids must match the
construction, and the Fortran Biot-Savart polygon of one group must equal
the closure that CoilsFile.write_simple applies to that group's filaments
(compared before its 15-digit text formatting).
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np


def load_oracle():
    path = Path(__file__).resolve().parents[2] / "python" / "libneo" / "coils.py"
    spec = importlib.util.spec_from_file_location("libneo_coils_oracle", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def ring(radius, z, phase, n=17):
    t = phase + 2.0 * np.pi * np.arange(n + 1) / n
    pts = np.column_stack([radius * np.cos(t), radius * np.sin(t),
                           z + 0.1 * np.sin(3 * t)])
    pts[-1] = pts[0]
    return pts


def tilted(phi, n=13):
    t = 2.0 * np.pi * np.arange(n + 1) / n
    r = 1.0 + 0.4 * np.cos(t)
    pts = np.column_stack([r * np.cos(phi + 0.05 * np.sin(t)),
                           r * np.sin(phi + 0.05 * np.sin(t)), 0.5 * np.sin(t)])
    pts[-1] = pts[0] + 1e-3
    return pts


def write_coil(lines, pts, current, group, name, repeat_first=False):
    body = pts[:-1]
    if repeat_first:
        body = np.vstack([body, body[:1]])
    for p in body:
        lines.append(f"{p[0]:.16e} {p[1]:.16e}\t{p[2]:.16e} {current:.16e}")
    term = f"{pts[-1][0]:.16e} {pts[-1][1]:.16e} {pts[-1][2]:.16e} 0.0"
    if group is not None:
        term += f" {group}" + (f" {name}" if name else "")
    lines.append(term + "  # terminator")


def synthetic(path: Path):
    lines = ["! synthetic MAKEGRID coils, nfp=3", "periods 3", "begin filament",
             "mirror NIL"]
    expected = []
    for k in range(3):
        write_coil(lines, tilted(2 * np.pi * k / 3), 108.0, 1, "NPC1")
        expected.append(1)
    write_coil(lines, ring(1.7, 0.6, 0.1), -36.0, 2, "PA", repeat_first=True)
    expected.append(2)
    lines.append("end filament")
    lines += ["", "# second block", "begin filament"]
    write_coil(lines, ring(1.7, -0.6, 0.2), -36.0, 2, None)
    expected.append(2)
    write_coil(lines, ring(0.8, 0.0, 0.0, n=9), 8.0, 4, "CC")
    expected.append(4)
    tail = ring(2.2, 0.0, 0.3, n=7)
    for p in tail[:-1]:
        lines.append(f"{p[0]:.16e} {p[1]:.16e} {p[2]:.16e} 5.0 ! inline")
    expected.append(0)
    lines.append("end")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return expected


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    periods = int(lines[0].split()[1])
    toks = lines[1].split()
    mirror = toks[1] if len(toks) > 1 else None
    n_coils, n_groups = map(int, lines[2].split()[1:3])
    coils, pos = [], 3
    for _ in range(n_coils):
        _, npts, group, name = lines[pos].split()
        data = np.array([[float(v) for v in ln.split()]
                         for ln in lines[pos + 1:pos + 1 + int(npts)]])
        coils.append((int(group), name, data))
        pos += 1 + int(npts)
    return periods, mirror, n_groups, coils


def read_point_list(path: Path):
    lines = path.read_text().split("\n")
    n = int(lines[0])
    return np.array([[float(v) for v in ln.split()] for ln in lines[1:1 + n]])


def fail(msg):
    print("FAIL:", msg)
    return 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()
    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    coils_path = work / "coils.makegrid_oracle"
    expected_groups = synthetic(coils_path)
    group = 2

    subprocess.run([args.exe, str(coils_path), str(work / "makegrid_dump.txt"),
                    str(work / "makegrid_group.txt"), str(group)], check=True)
    periods, mirror, n_groups, coils = read_dump(work / "makegrid_dump.txt")

    oracle = load_oracle()
    cf = oracle.CoilsFile.from_file(str(coils_path))
    if periods != cf.periods or mirror != cf.mirror:
        return fail(f"metadata {periods} {mirror} vs {cf.periods} {cf.mirror}")
    if len(coils) != len(cf.filaments):
        return fail(f"coil count {len(coils)} vs {len(cf.filaments)}")
    for k, ((g, _, data), fil) in enumerate(zip(coils, cf.filaments)):
        if not np.array_equal(data[:, :3], fil.coords):
            return fail(f"coil {k} coordinates differ")
        if not np.array_equal(data[:, 3], fil.current):
            return fail(f"coil {k} currents differ")
    groups = [g for g, _, _ in coils]
    if groups != expected_groups or n_groups != max(expected_groups):
        return fail(f"groups {groups} (n={n_groups}) vs {expected_groups}")
    names = {g: n for g, n, _ in coils}
    if names[1] != "NPC1" or names[2] != "PA" or names[4] != "CC":
        return fail(f"group names {names}")

    selected = [f for f, g in zip(cf.filaments, expected_groups) if g == group]
    closed = [oracle._closed_with_zero_end(f) for f in selected]
    py_pts = np.vstack([np.column_stack([c, i]) for c, i in closed])
    f_pts = read_point_list(work / "makegrid_group.txt")
    if py_pts.shape != f_pts.shape or not np.array_equal(py_pts, f_pts):
        return fail("group polygon differs from the write_simple closure")

    print(f"OK: {len(coils)} coils, groups {groups}, group {group} polygon "
          f"{f_pts.shape[0]} points identical to the Python oracle")
    return 0


if __name__ == "__main__":
    sys.exit(main())
