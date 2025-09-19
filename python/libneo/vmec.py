"""
VMEC geometry helpers: evaluate cylindrical coordinates (R, Z, phi)
from VMEC NetCDF outputs for given (s, theta, zeta).

Minimal API to support tests and downstream use. Focuses on read-only
coordinate evaluation using Fourier coefficients.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple
import numpy as np
from netCDF4 import Dataset


def _to_mode_ns_from_var(var) -> np.ndarray:
    """Ensure coefficient array layout is (nmode, ns) using NetCDF var dims.

    VMEC typically stores coefficients with dims ('radius','mn_mode') which is
    (ns, nmode). We transpose to (nmode, ns). If dims are already
    ('mn_mode','radius'), we leave as-is.
    """
    arr = np.array(var[:])
    if arr.ndim != 2:
        raise ValueError("Expected 2D coefficient array")
    dims = getattr(var, "dimensions", ())
    if len(dims) == 2 and dims[0] == "mn_mode" and dims[1] == "radius":
        return arr
    if len(dims) == 2 and dims[0] == "radius" and dims[1] == "mn_mode":
        return arr.T
    n0, n1 = arr.shape
    return arr.T if n0 < n1 else arr


def _cfunct(theta: np.ndarray, zeta: float, coeff: np.ndarray, xm: np.ndarray, xn: np.ndarray) -> np.ndarray:
    """Cosine series evaluation: sum_k coeff_k(s) * cos(xm_k*theta + xn_k*zeta)."""
    theta = np.asarray(theta, dtype=float)
    angle = np.outer(xm, theta) + np.outer(xn, np.atleast_1d(float(zeta)))
    cos_terms = np.cos(angle)
    return coeff.T @ cos_terms


def _sfunct(theta: np.ndarray, zeta: float, coeff: np.ndarray, xm: np.ndarray, xn: np.ndarray) -> np.ndarray:
    """Sine series evaluation: sum_k coeff_k(s) * sin(xm_k*theta + xn_k*zeta)."""
    theta = np.asarray(theta, dtype=float)
    angle = np.outer(xm, theta) + np.outer(xn, np.atleast_1d(float(zeta)))
    sin_terms = np.sin(angle)
    return coeff.T @ sin_terms


@dataclass
class VMECGeometry:
    xm: np.ndarray
    xn: np.ndarray
    rmnc: np.ndarray  # (nmode, ns)
    zmns: np.ndarray  # (nmode, ns)
    rmns: np.ndarray | None = None  # optional asymmetry
    zmnc: np.ndarray | None = None

    @classmethod
    def from_file(cls, nc_path: str) -> "VMECGeometry":
        with Dataset(nc_path, mode="r") as ds:
            xm = np.array(ds.variables["xm"][:])
            xn = np.array(ds.variables["xn"][:])
            rmnc = _to_mode_ns_from_var(ds.variables["rmnc"])
            zmns = _to_mode_ns_from_var(ds.variables["zmns"])

            rmns = zmnc = None
            lasym = False
            if "lasym__logical__" in ds.variables:
                lasym = bool(np.array(ds.variables["lasym__logical__"][...]))
            elif "lasym" in ds.variables:
                lasym = bool(np.array(ds.variables["lasym"][...]))
            if lasym:
                if "rmns" in ds.variables and "zmnc" in ds.variables:
                    rmns = _to_mode_ns_from_var(ds.variables["rmns"])
                    zmnc = _to_mode_ns_from_var(ds.variables["zmnc"])
        return cls(xm=xm, xn=xn, rmnc=rmnc, zmns=zmns, rmns=rmns, zmnc=zmnc)

    def coords(self, s_index: int, theta: np.ndarray, zeta: float, use_asym: bool = True) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Evaluate cylindrical coordinates (R, Z, phi) on a given s surface index.

        - s_index: integer surface index in [0, ns-1]
        - theta: array of poloidal angles (radians)
        - zeta: toroidal/geometric angle (radians)
        - use_asym: include asymmetric terms if available
        """
        R = _cfunct(theta, zeta, self.rmnc, self.xm, self.xn)[s_index, :]
        Z = _sfunct(theta, zeta, self.zmns, self.xm, self.xn)[s_index, :]
        if use_asym and self.rmns is not None and self.zmnc is not None:
            R = R + _sfunct(theta, zeta, self.rmns, self.xm, self.xn)[s_index, :]
            Z = Z + _cfunct(theta, zeta, self.zmnc, self.xm, self.xn)[s_index, :]
        return R, Z, float(zeta)


def vmec_to_cylindrical(nc_path: str, s_index: int, theta: np.ndarray, zeta: float, use_asym: bool = True) -> Tuple[np.ndarray, np.ndarray, float]:
    geom = VMECGeometry.from_file(nc_path)
    return geom.coords(s_index, theta, zeta, use_asym=use_asym)
