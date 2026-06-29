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

    @property
    def ns(self) -> int:
        """Number of radial surfaces in VMEC output."""
        return int(self.rmnc.shape[1])

    def _s_grid(self) -> np.ndarray:
        """Return the normalized toroidal flux grid s in [0, 1]."""
        ns = self.ns
        if self.phi is not None and self.phi.size == ns and float(self.phi[-1]) > 0.0:
            return np.asarray(self.phi, dtype=float) / float(self.phi[-1])
        return np.linspace(0.0, 1.0, ns, dtype=float)

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
            if denom == 0.0:
                alpha = 0.0
            else:
                num = float(phi_target - self.phi[i0])
                if num == 0.0 and phi_target > float(self.phi[i0]):
                    num = float(np.nextafter(0.0, 1.0))
                alpha = float(num / denom)
                alpha = min(1.0, max(0.0, alpha))
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

    def coords_s_with_deriv(
        self,
        s: float,
        theta: np.ndarray,
        zeta: float,
        use_asym: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate (R, Z) and their radial derivatives (dR/ds, dZ/ds) at surface s.

        Uses second-order finite differences on the VMEC radial grid.
        At s=1, uses backward difference; at s=0, uses forward difference.

        Returns: (R, Z, dR_ds, dZ_ds) where each is an array over theta.
        """
        s_val = float(s)
        if not (0.0 <= s_val <= 1.0):
            raise ValueError("s must be in [0, 1]")

        R, Z, _ = self.coords_s(s_val, theta, zeta, use_asym=use_asym)

        s_grid = self._s_grid()
        ns = self.ns

        if ns < 3:
            dR_ds = np.zeros_like(R)
            dZ_ds = np.zeros_like(Z)
            return R, Z, dR_ds, dZ_ds

        if s_val >= 1.0:
            s0, s1, s2 = float(s_grid[-3]), float(s_grid[-2]), float(s_grid[-1])
            R0, Z0, _ = self.coords_s(s0, theta, zeta, use_asym=use_asym)
            R1, Z1, _ = self.coords_s(s1, theta, zeta, use_asym=use_asym)
            R2, Z2, _ = self.coords_s(s2, theta, zeta, use_asym=use_asym)
            h1, h2 = s1 - s0, s2 - s1
            dR_ds = (h1**2 * R2 - (h1**2 - h2**2) * R1 - h2**2 * R0) / (h1 * h2 * (h1 + h2))
            dZ_ds = (h1**2 * Z2 - (h1**2 - h2**2) * Z1 - h2**2 * Z0) / (h1 * h2 * (h1 + h2))
        elif s_val <= 0.0:
            s0, s1, s2 = float(s_grid[0]), float(s_grid[1]), float(s_grid[2])
            R0, Z0, _ = self.coords_s(s0, theta, zeta, use_asym=use_asym)
            R1, Z1, _ = self.coords_s(s1, theta, zeta, use_asym=use_asym)
            R2, Z2, _ = self.coords_s(s2, theta, zeta, use_asym=use_asym)
            h1, h2 = s1 - s0, s2 - s1
            dR_ds = (-(h1 + h2)**2 * R0 + (h1 + h2)**2 * R1 - h1**2 * R1 + h1**2 * R2) / (h1 * h2 * (h1 + h2))
            dZ_ds = (-(h1 + h2)**2 * Z0 + (h1 + h2)**2 * Z1 - h1**2 * Z1 + h1**2 * Z2) / (h1 * h2 * (h1 + h2))
        else:
            idx = int(np.searchsorted(s_grid, s_val, side="right"))
            idx = max(1, min(ns - 2, idx))
            s0, s1, s2 = float(s_grid[idx - 1]), float(s_grid[idx]), float(s_grid[idx + 1])
            R0, Z0, _ = self.coords_s(s0, theta, zeta, use_asym=use_asym)
            R1, Z1, _ = self.coords_s(s1, theta, zeta, use_asym=use_asym)
            R2, Z2, _ = self.coords_s(s2, theta, zeta, use_asym=use_asym)
            w0 = (2.0 * s_val - s1 - s2) / ((s0 - s1) * (s0 - s2))
            w1 = (2.0 * s_val - s0 - s2) / ((s1 - s0) * (s1 - s2))
            w2 = (2.0 * s_val - s0 - s1) / ((s2 - s0) * (s2 - s1))
            dR_ds = w0 * R0 + w1 * R1 + w2 * R2
            dZ_ds = w0 * Z0 + w1 * Z1 + w2 * Z2

        return R, Z, dR_ds, dZ_ds

    def coords_s_with_second_deriv(
        self,
        s: float,
        theta: np.ndarray,
        zeta: float,
        use_asym: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate (R, Z), first and second radial derivatives at surface s.

        Returns: (R, Z, dR_ds, dZ_ds, d2R_ds2, d2Z_ds2) where each is array over theta.
        """
        s_val = float(s)
        if not (0.0 <= s_val <= 1.0):
            raise ValueError("s must be in [0, 1]")

        R, Z, dR_ds, dZ_ds = self.coords_s_with_deriv(s_val, theta, zeta, use_asym)

        s_grid = self._s_grid()
        ns = self.ns

        if ns < 4:
            d2R_ds2 = np.zeros_like(R)
            d2Z_ds2 = np.zeros_like(Z)
            return R, Z, dR_ds, dZ_ds, d2R_ds2, d2Z_ds2

        if s_val >= 1.0:
            s0, s1, s2 = float(s_grid[-3]), float(s_grid[-2]), float(s_grid[-1])
            _, _, dR0, dZ0 = self.coords_s_with_deriv(s0, theta, zeta, use_asym)
            _, _, dR1, dZ1 = self.coords_s_with_deriv(s1, theta, zeta, use_asym)
            _, _, dR2, dZ2 = self.coords_s_with_deriv(s2, theta, zeta, use_asym)
            h1, h2 = s1 - s0, s2 - s1
            d2R_ds2 = (h1**2 * dR2 - (h1**2 - h2**2) * dR1 - h2**2 * dR0) / (
                h1 * h2 * (h1 + h2)
            )
            d2Z_ds2 = (h1**2 * dZ2 - (h1**2 - h2**2) * dZ1 - h2**2 * dZ0) / (
                h1 * h2 * (h1 + h2)
            )
        elif s_val <= 0.0:
            s0, s1, s2 = float(s_grid[0]), float(s_grid[1]), float(s_grid[2])
            _, _, dR0, dZ0 = self.coords_s_with_deriv(s0, theta, zeta, use_asym)
            _, _, dR1, dZ1 = self.coords_s_with_deriv(s1, theta, zeta, use_asym)
            _, _, dR2, dZ2 = self.coords_s_with_deriv(s2, theta, zeta, use_asym)
            h1, h2 = s1 - s0, s2 - s1
            d2R_ds2 = (
                -(h1 + h2) ** 2 * dR0
                + (h1 + h2) ** 2 * dR1
                - h1**2 * dR1
                + h1**2 * dR2
            ) / (h1 * h2 * (h1 + h2))
            d2Z_ds2 = (
                -(h1 + h2) ** 2 * dZ0
                + (h1 + h2) ** 2 * dZ1
                - h1**2 * dZ1
                + h1**2 * dZ2
            ) / (h1 * h2 * (h1 + h2))
        else:
            idx = int(np.searchsorted(s_grid, s_val, side="right"))
            idx = max(1, min(ns - 2, idx))
            s0, s1, s2 = float(s_grid[idx - 1]), float(s_grid[idx]), float(s_grid[idx + 1])
            _, _, dR0, dZ0 = self.coords_s_with_deriv(s0, theta, zeta, use_asym)
            _, _, dR1, dZ1 = self.coords_s_with_deriv(s1, theta, zeta, use_asym)
            _, _, dR2, dZ2 = self.coords_s_with_deriv(s2, theta, zeta, use_asym)
            w0 = (2.0 * s_val - s1 - s2) / ((s0 - s1) * (s0 - s2))
            w1 = (2.0 * s_val - s0 - s2) / ((s1 - s0) * (s1 - s2))
            w2 = (2.0 * s_val - s0 - s1) / ((s2 - s0) * (s2 - s1))
            d2R_ds2 = w0 * dR0 + w1 * dR1 + w2 * dR2
            d2Z_ds2 = w0 * dZ0 + w1 * dZ1 + w2 * dZ2

        return R, Z, dR_ds, dZ_ds, d2R_ds2, d2Z_ds2

    def boundary_rz(
        self,
        s: float,
        theta: np.ndarray,
        zeta: float,
        *,
        boundary_offset: float = 0.0,
        use_asym: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Evaluate a VMEC (R,Z) boundary curve with an optional normal offset.

        - `boundary_offset` adds a fixed outward offset in meters along the local
          outward curve normal in the R-Z plane at fixed zeta.
        """
        boundary_offset = float(boundary_offset)
        if boundary_offset < 0.0:
            raise ValueError("boundary_offset must be >= 0")

        if boundary_offset == 0.0:
            R, Z, _ = self.coords_s(float(s), theta, float(zeta), use_asym=use_asym)
            return R, Z, float(zeta)

        try:
            from shapely.geometry import Polygon
        except Exception as exc:  # pragma: no cover
            raise ImportError(
                "boundary_offset requires shapely; install libneo with the 'chartmap' extra"
            ) from exc

        th = np.asarray(theta, dtype=float)
        if th.ndim == 0:
            th = th.reshape((1,))
        if th.ndim != 1:
            raise ValueError("theta must be a scalar or 1D array")

        n_base = 4096
        th_base = np.linspace(0.0, 2.0 * np.pi, n_base, endpoint=False, dtype=float)
        R_base, Z_base, _ = self.coords_s(float(s), th_base, float(zeta), use_asym=use_asym)
        poly = Polygon(np.column_stack([R_base, Z_base]))
        if not poly.is_valid:
            poly = poly.buffer(0.0)
        if poly.is_empty:
            raise ValueError("invalid VMEC boundary polygon for buffering")

        off = poly.buffer(boundary_offset, join_style=1)
        if off.is_empty:
            raise ValueError("boundary_offset produced empty geometry")
        if off.geom_type != "Polygon":
            off = max(list(off.geoms), key=lambda g: g.area)

        coords = np.asarray(off.exterior.coords, dtype=float)
        if coords.ndim != 2 or coords.shape[0] < 4:
            raise ValueError("unexpected buffered boundary geometry")

        # Canonicalize startpoint and orientation so theta=0 is deterministic.
        coords_open = coords[:-1, :]
        i0 = int(np.argmax(coords_open[:, 0]))
        coords_open = np.vstack([coords_open[i0:, :], coords_open[:i0, :]])
        if coords_open.shape[0] >= 2 and coords_open[1, 1] - coords_open[0, 1] < 0.0:
            coords_open = coords_open[::-1, :]
        coords = np.vstack([coords_open, coords_open[0, :]])

        seg = np.sqrt(np.sum((coords[1:] - coords[:-1]) ** 2, axis=1))
        s_coords = np.concatenate(([0.0], np.cumsum(seg)))
        length = float(s_coords[-1])
        if length == 0.0:
            raise ValueError("buffered boundary has zero length")

        u = (th % (2.0 * np.pi)) / (2.0 * np.pi)
        s_query = u * length
        R_out = np.interp(s_query, s_coords, coords[:, 0])
        Z_out = np.interp(s_query, s_coords, coords[:, 1])
        return R_out, Z_out, float(zeta)

    def find_theta(
        self,
        R_target: float | np.ndarray,
        Z_target: float | np.ndarray,
        zeta: float | np.ndarray,
        s: float = 1.0,
        *,
        theta_guess: float | np.ndarray | None = None,
        max_iter: int = 20,
        tol: float = 1.0e-10,
        n_starts: int = 12,
        use_asym: bool = True,
    ) -> np.ndarray:
        """
        Find VMEC poloidal angle theta for given (R, Z) at fixed (s, zeta).

        Uses Newton-Raphson iteration with analytic Jacobian from Fourier
        derivatives. Falls back to multi-start if no initial guess is provided.

        Parameters
        ----------
        R_target : target cylindrical R coordinate(s) in meters
        Z_target : target cylindrical Z coordinate(s) in meters
        zeta : toroidal angle(s) in radians
        s : normalized toroidal flux coordinate (default 1.0 = LCFS)
        theta_guess : initial guess for theta (radians), or None for multi-start
        max_iter : maximum Newton iterations per start
        tol : convergence tolerance (squared error in R, Z)
        n_starts : number of initial guesses for multi-start (if no guess given)
        use_asym : include asymmetric Fourier terms if available

        Returns
        -------
        theta : VMEC poloidal angle(s) in [0, 2*pi)
        """
        R_arr = np.atleast_1d(np.asarray(R_target, dtype=float)).flatten()
        Z_arr = np.atleast_1d(np.asarray(Z_target, dtype=float)).flatten()
        zeta_arr = np.atleast_1d(np.asarray(zeta, dtype=float)).flatten()

        if not (R_arr.shape == Z_arr.shape == zeta_arr.shape):
            raise ValueError("R_target, Z_target, zeta must have the same shape")

        n = R_arr.size
        theta_out = np.zeros(n, dtype=float)
        s_val = float(s)
        TWOPI = 2.0 * np.pi

        # Precompute coefficients at surface s for efficiency
        ns = self.ns
        if ns == 1:
            rmnc_s = self.rmnc[:, 0]
            zmns_s = self.zmns[:, 0]
            rmns_s = self.rmns[:, 0] if self.rmns is not None else None
            zmnc_s = self.zmnc[:, 0] if self.zmnc is not None else None
        else:
            x = s_val * float(ns - 1)
            i0 = int(np.floor(x))
            i0 = max(0, min(ns - 2, i0))
            i1 = i0 + 1
            alpha = x - float(i0)
            rmnc_s = (1.0 - alpha) * self.rmnc[:, i0] + alpha * self.rmnc[:, i1]
            zmns_s = (1.0 - alpha) * self.zmns[:, i0] + alpha * self.zmns[:, i1]
            if use_asym and self.rmns is not None and self.zmnc is not None:
                rmns_s = (1.0 - alpha) * self.rmns[:, i0] + alpha * self.rmns[:, i1]
                zmnc_s = (1.0 - alpha) * self.zmnc[:, i0] + alpha * self.zmnc[:, i1]
            else:
                rmns_s = None
                zmnc_s = None

        def _eval_rz(th: float, z: float) -> Tuple[float, float, float, float]:
            """Evaluate (R, Z, dR/dtheta, dZ/dtheta) at single point."""
            angle = self.xm * th - self.xn * z
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            R = float(np.dot(rmnc_s, cos_a))
            Z = float(np.dot(zmns_s, sin_a))
            dR_dth = float(np.dot(-self.xm * rmnc_s, sin_a))
            dZ_dth = float(np.dot(self.xm * zmns_s, cos_a))
            if rmns_s is not None and zmnc_s is not None:
                R += float(np.dot(rmns_s, sin_a))
                Z += float(np.dot(zmnc_s, cos_a))
                dR_dth += float(np.dot(self.xm * rmns_s, cos_a))
                dZ_dth += float(np.dot(-self.xm * zmnc_s, sin_a))
            return R, Z, dR_dth, dZ_dth

        def _newton(th0: float, z: float, R_t: float, Z_t: float) -> Tuple[float, float]:
            """Newton iteration from initial guess, returns (theta, error)."""
            th = th0
            for _ in range(max_iter):
                R, Z, dR, dZ = _eval_rz(th, z)
                err_R = R - R_t
                err_Z = Z - Z_t
                err2 = err_R * err_R + err_Z * err_Z
                if err2 < tol:
                    return th % TWOPI, err2
                # Gauss-Newton step: minimize ||f||^2
                grad = 2.0 * (err_R * dR + err_Z * dZ)
                hess = 2.0 * (dR * dR + dZ * dZ)
                if abs(hess) < 1.0e-14:
                    break
                th -= grad / hess
            R, Z, _, _ = _eval_rz(th, z)
            err2 = (R - R_t) ** 2 + (Z - Z_t) ** 2
            return th % TWOPI, err2

        # Process each point
        if theta_guess is not None:
            guess_arr = np.atleast_1d(np.asarray(theta_guess, dtype=float)).flatten()
            if guess_arr.size == 1:
                guess_arr = np.full(n, guess_arr[0], dtype=float)
            elif guess_arr.size != n:
                raise ValueError("theta_guess must be scalar or match input size")
        else:
            guess_arr = None

        for i in range(n):
            R_t, Z_t, z = float(R_arr[i]), float(Z_arr[i]), float(zeta_arr[i])
            if guess_arr is not None:
                # Single Newton from provided guess
                th, _ = _newton(float(guess_arr[i]), z, R_t, Z_t)
                theta_out[i] = th
            else:
                # Multi-start: try n_starts initial guesses
                best_th = 0.0
                best_err = float("inf")
                for k in range(n_starts):
                    th0 = TWOPI * k / n_starts
                    th, err2 = _newton(th0, z, R_t, Z_t)
                    if err2 < best_err:
                        best_err = err2
                        best_th = th
                theta_out[i] = best_th

        return theta_out


def vmec_to_cylindrical(nc_path: str, s_index: int, theta: np.ndarray, zeta: float, use_asym: bool = True) -> Tuple[np.ndarray, np.ndarray, float]:
    geom = VMECGeometry.from_file(nc_path)
    return geom.coords(s_index, theta, zeta, use_asym=use_asym)
