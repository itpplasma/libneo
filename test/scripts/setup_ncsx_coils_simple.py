#!/usr/bin/env python3
"""
Download an NCSX coils file and convert it to SIMPLE coils format.

Outputs in the specified directory:
- ncsx_coils.c09r00
- ncsx_coils.simple
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys
import urllib.request


COILS_URL = "https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00"
COILS_SHA256 = "3c429da06f4c062887a497a16e2d2bd10f0ecb0b8858c252698631f3853da428"


def _import_libneo_from_repo() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "python"))


def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(path: Path) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    try:
        with urllib.request.urlopen(COILS_URL, timeout=120) as response, tmp_path.open(
            "wb"
        ) as handle:
            handle.write(response.read())
        got = sha256sum(tmp_path)
        if got != COILS_SHA256:
            raise RuntimeError(f"sha256 mismatch: expected {COILS_SHA256}, got {got}")
        tmp_path.replace(path)
    finally:
        tmp_path.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_ncsx_coils_simple")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    coils_src = out_dir / "ncsx_coils.c09r00"
    coils_simple = out_dir / "ncsx_coils.simple"

    if coils_src.exists() and coils_src.stat().st_size > 0:
        got = sha256sum(coils_src)
        if got != COILS_SHA256:
            coils_src.unlink()
        else:
            _import_libneo_from_repo()
            from libneo.coils import CoilsFile

            cf = CoilsFile.from_file(str(coils_src))
            cf.write_simple(str(coils_simple))
            if not coils_simple.exists() or coils_simple.stat().st_size == 0:
                raise RuntimeError("failed to write SIMPLE coils file")
            return 0

    download(coils_src)
    _import_libneo_from_repo()
    from libneo.coils import CoilsFile

    cf = CoilsFile.from_file(str(coils_src))
    cf.write_simple(str(coils_simple))
    if not coils_simple.exists() or coils_simple.stat().st_size == 0:
        raise RuntimeError("failed to write SIMPLE coils file")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

