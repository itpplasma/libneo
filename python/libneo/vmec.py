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
    """Cosine series evaluation: sum_k coeff_k(s) * cos(xm_k*theta - xn_k*zeta)."""
    theta = np.asarray(theta, dtype=float)
    angle = np.outer(xm, theta) - np.outer(xn, np.atleast_1d(float(zeta)))
    cos_terms = np.cos(angle)
    return coeff.T @ cos_terms


def _sfunct(theta: np.ndarray, zeta: float, coeff: np.ndarray, xm: np.ndarray, xn: np.ndarray) -> np.ndarray:
    """Sine series evaluation: sum_k coeff_k(s) * sin(xm_k*theta - xn_k*zeta)."""
    theta = np.asarray(theta, dtype=float)
    angle = np.outer(xm, theta) - np.outer(xn, np.atleast_1d(float(zeta)))
    sin_terms = np.sin(angle)
    return coeff.T @ sin_terms


@dataclass
class VMECGeometry:
    xm: np.ndarray
    xn: np.ndarray
    phi: np.ndarray | None  # (ns,) toroidal flux coordinate, if present
    rmnc: np.ndarray  # (nmode, ns)
    zmns: np.ndarray  # (nmode, ns)
    rmns: np.ndarray | None = None  # optional asymmetry
    zmnc: np.ndarray | None = None

    @classmethod
    def from_file(cls, nc_path: str) -> "VMECGeometry":
        with Dataset(nc_path, mode="r") as ds:
            xm = np.array(ds.variables["xm"][:])
            xn = np.array(ds.variables["xn"][:])
            phi = None
            if "phi" in ds.variables:
                phi = np.array(ds.variables["phi"][:], dtype=float)
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
        return cls(xm=xm, xn=xn, phi=phi, rmnc=rmnc, zmns=zmns, rmns=rmns, zmnc=zmnc)

    def coords(self, s_index: int, theta: np.ndarray, zeta: float, use_asym: bool = True) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Evaluate cylindrical coordinates (R, Z, phi) on a given s surface index.

        - s_index: integer surface index in [0, ns-1]
        - theta: array of poloidal angles (radians)
        - zeta: toroidal/geometric angle (radians, full-torus cylindrical angle)
        - use_asym: include asymmetric terms if available

        Notes
        -----
        VMEC `wout_*.nc` files store `xn = n*nfp` (signed), so the Fourier phase is
        evaluated as `m*theta - xn*zeta`. No additional `nfp` factor should be
        applied to `zeta` here.
        """
        R = _cfunct(theta, zeta, self.rmnc, self.xm, self.xn)[s_index, :]
        Z = _sfunct(theta, zeta, self.zmns, self.xm, self.xn)[s_index, :]
        if use_asym and self.rmns is not None and self.zmnc is not None:
            R = R + _sfunct(theta, zeta, self.rmns, self.xm, self.xn)[s_index, :]
            Z = Z + _cfunct(theta, zeta, self.zmnc, self.xm, self.xn)[s_index, :]
        return R, Z, float(zeta)

    def coords_s(self, s: float, theta: np.ndarray, zeta: float, use_asym: bool = True) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Evaluate cylindrical coordinates (R, Z, phi) on a fractional surface coordinate s in [0, 1].

        If the wout file provides the `phi` radial coordinate array, s is interpreted as
        normalized toroidal flux: s = phi / phi_edge. Otherwise, s is mapped linearly
        to the discrete surface index range [0, ns-1].
        """
        s_val = float(s)
        if not (0.0 <= s_val <= 1.0):
            raise ValueError("s must be in [0, 1]")

        ns = int(self.rmnc.shape[1])
        if ns < 1:
            raise ValueError("invalid VMEC coefficient arrays: ns must be >= 1")

        if ns == 1 or s_val == 0.0:
            return self.coords(0, theta, zeta, use_asym=use_asym)
        if s_val == 1.0:
            return self.coords(ns - 1, theta, zeta, use_asym=use_asym)

        if self.phi is not None and self.phi.size == ns and float(self.phi[-1]) > 0.0:
            phi_edge = float(self.phi[-1])
            phi_target = s_val * phi_edge
            i1 = int(np.searchsorted(self.phi, phi_target, side="right"))
            i1 = max(1, min(ns - 1, i1))
            i0 = i1 - 1
            denom = float(self.phi[i1] - self.phi[i0])
            alpha = 0.0 if denom == 0.0 else float((phi_target - self.phi[i0]) / denom)
        else:
            x = s_val * float(ns - 1)
            i0 = int(np.floor(x))
            i0 = max(0, min(ns - 2, i0))
            i1 = i0 + 1
            alpha = float(x - float(i0))

        w0 = 1.0 - alpha
        w1 = alpha
        rmnc_s = (w0 * self.rmnc[:, i0] + w1 * self.rmnc[:, i1])[:, None]
        zmns_s = (w0 * self.zmns[:, i0] + w1 * self.zmns[:, i1])[:, None]

        R = _cfunct(theta, zeta, rmnc_s, self.xm, self.xn)[0, :]
        Z = _sfunct(theta, zeta, zmns_s, self.xm, self.xn)[0, :]
        if use_asym and self.rmns is not None and self.zmnc is not None:
            rmns_s = (w0 * self.rmns[:, i0] + w1 * self.rmns[:, i1])[:, None]
            zmnc_s = (w0 * self.zmnc[:, i0] + w1 * self.zmnc[:, i1])[:, None]
            R = R + _sfunct(theta, zeta, rmns_s, self.xm, self.xn)[0, :]
            Z = Z + _cfunct(theta, zeta, zmnc_s, self.xm, self.xn)[0, :]
        return R, Z, float(zeta)


def vmec_to_cylindrical(nc_path: str, s_index: int, theta: np.ndarray, zeta: float, use_asym: bool = True) -> Tuple[np.ndarray, np.ndarray, float]:
    geom = VMECGeometry.from_file(nc_path)
    return geom.coords(s_index, theta, zeta, use_asym=use_asym)
