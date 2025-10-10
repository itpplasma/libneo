#!/usr/bin/env python3
"""Download a public GEQDSK sample and invoke test_geoflux.x."""

from __future__ import annotations

import argparse
import os
import pathlib
import shutil
import ssl
import subprocess
import sys
import urllib.request

GEQDSK_URL = "https://crppwww.epfl.ch/~sauter/benchmark/EQDSK_I"
SSL_CONTEXT = ssl._create_unverified_context()
LIBNEO_TESTING_ENV = "LIBNEO_TESTING"
_DISABLED_VALUES = {"", "0", "false", "no", "off"}


def testing_enabled() -> bool:
    value = os.environ.get(LIBNEO_TESTING_ENV, "")
    return value.strip().lower() not in _DISABLED_VALUES


def require_testing_enabled(action: str) -> None:
    if testing_enabled():
        return
    print(
        f"{action} skipped: set {LIBNEO_TESTING_ENV}=1 to enable", file=sys.stderr
    )
    raise SystemExit(0)


def download(url: str, destination: pathlib.Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = destination.with_suffix(destination.suffix + ".tmp")
    try:
        with urllib.request.urlopen(url, context=SSL_CONTEXT) as response, \
                tmp_path.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        tmp_path.replace(destination)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def ensure_geqdsk(dest_dir: pathlib.Path) -> pathlib.Path:
    target = dest_dir / "EQDSK_I.geqdsk"
    if not target.exists() or target.stat().st_size == 0:
        download(GEQDSK_URL, target)
    return target


def run_test(executable: pathlib.Path, geqdsk_path: pathlib.Path) -> None:
    result = subprocess.run([str(executable), str(geqdsk_path)], check=False)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, result.args, result.returncode)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="run_geoflux_test")
    parser.add_argument("--exe", required=True, help="Path to test_geoflux.x")
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory to store downloaded GEQDSK files",
    )
    args = parser.parse_args(argv)

    require_testing_enabled("geoflux integration test")

    exe_path = pathlib.Path(args.exe).resolve()
    if not exe_path.is_file():
        parser.error(f"executable not found: {exe_path}")

    data_dir = pathlib.Path(args.data_dir).resolve()

    geqdsk_path = ensure_geqdsk(data_dir)
    run_test(exe_path, geqdsk_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
