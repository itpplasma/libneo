"""Run hiddenSymmetries booz_xform on a pinned VMEC wout file.

The wrapper makes the Fourier resolution and selected half-grid surfaces
explicit and records enough provenance in the resulting boozmn NetCDF file to
identify the exact VMEC input used for the transform.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _record_provenance(
    output: Path,
    source: Path,
    *,
    version: str,
    mboz: int,
    nboz: int,
    surfaces: np.ndarray,
) -> None:
    import netCDF4

    with netCDF4.Dataset(source) as vmec:
        aminor_m = float(np.asarray(vmec.variables["Aminor_p"][:]))
        rmajor_m = float(np.asarray(vmec.variables["Rmajor_p"][:]))

    with netCDF4.Dataset(output, "a") as dataset:
        dataset.libneo_booz_xform_source = str(source.resolve())
        dataset.libneo_booz_xform_source_sha256 = _sha256(source)
        dataset.libneo_booz_xform_version = version
        dataset.libneo_booz_xform_mboz = mboz
        dataset.libneo_booz_xform_nboz = nboz
        dataset.libneo_booz_xform_surfaces = ",".join(
            str(surface) for surface in surfaces
        )
        dataset.libneo_vmec_aminor_m = aminor_m
        dataset.libneo_vmec_rmajor_m = rmajor_m


def run_booz_xform(
    wout,
    output,
    *,
    mboz=64,
    nboz=0,
    surfaces=None,
):
    """Transform a VMEC equilibrium and return the output path.

    ``surfaces`` contains zero-based indices into booz_xform's VMEC half grid.
    Omitting it transforms every available half-grid surface.
    """
    import booz_xform

    wout = Path(wout)
    output = Path(output)
    if mboz < 1:
        raise ValueError("mboz must be positive")
    if nboz < 0:
        raise ValueError("nboz must be nonnegative")

    transform = booz_xform.Booz_xform()
    # flux=True is essential: without it booz_xform writes phi_b=0 and the
    # resulting file cannot preserve the VMEC toroidal-flux normalization.
    transform.read_wout(str(wout), flux=True)
    transform.mboz = int(mboz)
    transform.nboz = int(nboz)

    if surfaces is None:
        selected = np.asarray(transform.compute_surfs, dtype=int)
    else:
        selected = np.asarray(surfaces, dtype=int)
        if selected.ndim != 1 or selected.size == 0:
            raise ValueError("surfaces must be a nonempty one-dimensional list")
        if np.any(selected < 0) or np.any(selected >= transform.ns_in):
            raise ValueError(
                f"surfaces outside zero-based half-grid range "
                f"0..{transform.ns_in - 1}"
            )
        if len(np.unique(selected)) != len(selected):
            raise ValueError("surfaces must not contain duplicates")
        transform.compute_surfs = selected

    transform.run()
    output.parent.mkdir(parents=True, exist_ok=True)
    transform.write_boozmn(str(output))
    _record_provenance(
        output,
        wout,
        version=str(getattr(booz_xform, "__version__", "unknown")),
        mboz=int(mboz),
        nboz=int(nboz),
        surfaces=selected,
    )
    return output


def _parse_surfaces(value):
    if value is None:
        return None
    try:
        return [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "surfaces must be comma-separated integers"
        ) from error


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run booz_xform on a VMEC wout file with provenance"
    )
    parser.add_argument("wout", help="VMEC wout NetCDF file")
    parser.add_argument("output", help="Output boozmn NetCDF file")
    parser.add_argument("--mboz", type=int, default=64)
    parser.add_argument("--nboz", type=int, default=0)
    parser.add_argument(
        "--surfaces",
        type=_parse_surfaces,
        help="comma-separated zero-based VMEC half-grid indices (default: all)",
    )
    args = parser.parse_args(argv)
    output = run_booz_xform(
        args.wout,
        args.output,
        mboz=args.mboz,
        nboz=args.nboz,
        surfaces=args.surfaces,
    )
    print(output)


if __name__ == "__main__":
    main()
