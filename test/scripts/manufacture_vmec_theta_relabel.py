#!/usr/bin/env python3
"""Create the pure VMEC coordinate relabel ``theta'=-theta``."""

from __future__ import annotations

import shutil
from pathlib import Path

import netCDF4


SCALAR_SINE = {"rmns", "zmns", "bmns"}
ANGLE_COSINE = {"lmnc"}
ORIENTED_DENSITY_COSINE = {"gmnc"}
THETA_COMPONENT_COSINE = {"bsubumnc", "bsupumnc", "currumnc"}
UNCHANGED_COMPONENT_SINE = {"bsubvmns", "bsupvmns", "currvmns", "bsubsmns"}


def negate_if_present(dataset: netCDF4.Dataset, names: set[str]) -> None:
    for name in names:
        if name in dataset.variables:
            dataset.variables[name][:] = -dataset.variables[name][:]


def relabel_wout(source: Path, output: Path) -> None:
    """Write the same Cartesian equilibrium in the ``theta'=-theta`` chart."""

    shutil.copyfile(source, output)
    with netCDF4.Dataset(output, "a") as dataset:
        signgs = int(dataset.variables["signgs"][:])
        if signgs not in (-1, 1):
            raise RuntimeError(f"invalid source signgs={signgs}")
        dataset.variables["signgs"].assignValue(-signgs)
        for name in ("xn", "xn_nyq"):
            dataset.variables[name][:] = -dataset.variables[name][:]
        for name in ("iotaf", "iotas", "q_factor", "chi", "chipf"):
            if name in dataset.variables:
                dataset.variables[name][:] = -dataset.variables[name][:]

        negate_if_present(dataset, SCALAR_SINE)
        negate_if_present(dataset, ANGLE_COSINE)
        negate_if_present(dataset, ORIENTED_DENSITY_COSINE)
        negate_if_present(dataset, THETA_COMPONENT_COSINE)
        negate_if_present(dataset, UNCHANGED_COMPONENT_SINE)
        for name in ("buco", "jcuru"):
            if name in dataset.variables:
                dataset.variables[name][:] = -dataset.variables[name][:]

        dataset.setncattr(
            "libneo_manufactured_relabel",
            "theta_prime=-theta; zeta_prime=zeta; pure coordinate relabel",
        )
