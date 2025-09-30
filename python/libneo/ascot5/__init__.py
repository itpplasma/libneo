"""ASCOT5 field conversion helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import math
import time

import h5py
import netCDF4
import numpy as np

from libneo.vmec import VMECGeometry

try:  # pragma: no cover - exercised in integration tests
    import _magfie  # type: ignore

    _VMEC_WRAPPERS = _magfie.f2py_vmec_wrappers
except (ImportError, AttributeError) as exc:  # pragma: no cover
    raise RuntimeError(
        "_magfie extension with VMEC wrappers is required. "
        "Ensure libneo is built with the Fortran interface." 
    ) from exc


CM_TO_M = 1.0e-2
GAUSS_TO_TESLA = 1.0e-4


@dataclass(frozen=True)
class B3DSField:
    """Container describing an ASCOT5 B_3DS field."""

    r_grid: np.ndarray
    z_grid: np.ndarray
    phi_grid_deg: np.ndarray
    axis_r: float
    axis_z: float
    psi: np.ndarray
    psi0: float
    psi1: float
    br: np.ndarray
    bphi: np.ndarray
    bz: np.ndarray

    def as_dict(self) -> Dict[str, object]:
        return {
            "b_rmin": float(self.r_grid[0]),
            "b_rmax": float(self.r_grid[-1]),
            "b_nr": int(self.r_grid.size),
            "b_zmin": float(self.z_grid[0]),
            "b_zmax": float(self.z_grid[-1]),
            "b_nz": int(self.z_grid.size),
            "b_phimin": float(self.phi_grid_deg[0]),
            "b_phimax": float(self.phi_grid_deg[-1] + (self.phi_grid_deg[1] - self.phi_grid_deg[0]) if self.phi_grid_deg.size > 1 else self.phi_grid_deg[0]),
            "b_nphi": int(self.phi_grid_deg.size),
            "axisr": float(self.axis_r),
            "axisz": float(self.axis_z),
            "psi": self.psi,
            "psi0": float(self.psi0),
            "psi1": float(self.psi1),
            "br": self.br,
            "bphi": self.bphi,
            "bz": self.bz,
        }


def _load_vmec_metadata(wout_path: Path) -> Tuple[int, float, float, np.ndarray]:
    with netCDF4.Dataset(wout_path) as ds:
        nfp = int(ds.variables["nfp"][()])
        axis_r = float(ds.variables["raxis_cc"][0])
        if "zaxis_cs" in ds.variables:
            axis_z = float(ds.variables["zaxis_cs"][0])
        else:
            axis_z = float(ds.variables["zaxis_cc"][0])
        iota = np.array(ds.variables["iotaf"][...], dtype=float)
    return nfp, axis_r, axis_z, iota


def _prepare_outer_surface(geom: VMECGeometry, n_phi: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta_samples = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
    phi_samples = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    ns_index = geom.rmnc.shape[1] - 1

    r_list = []
    z_list = []
    for phi in phi_samples:
        r_vals, z_vals, _ = geom.coords(ns_index, theta_samples, phi, use_asym=True)
        r_list.append(r_vals)
        z_list.append(z_vals)
    r_arr = np.stack(r_list)  # (nphi, ntheta)
    z_arr = np.stack(z_list)
    return phi_samples, r_arr, z_arr


def _initial_guess(
    axis_r: float,
    axis_z: float,
    phi_idx: int,
    phi_samples: np.ndarray,
    r_edge: np.ndarray,
    z_edge: np.ndarray,
    r_target: float,
    z_target: float,
) -> Tuple[float, float]:
    theta_guess = math.atan2(z_target - axis_z, r_target - axis_r)
    theta_guess = (theta_guess + 2.0 * math.pi) % (2.0 * math.pi)

    theta_samples = np.linspace(0.0, 2.0 * np.pi, r_edge.shape[1], endpoint=False)
    theta_ext = np.concatenate([theta_samples, theta_samples + 2.0 * np.pi])
    r_edge_ext = np.concatenate([r_edge[phi_idx], r_edge[phi_idx]])
    z_edge_ext = np.concatenate([z_edge[phi_idx], z_edge[phi_idx]])
    r_edge_val = np.interp(theta_guess, theta_ext, r_edge_ext)
    z_edge_val = np.interp(theta_guess, theta_ext, z_edge_ext)

    axis_vec = np.array([axis_r, axis_z])
    target_vec = np.array([r_target, z_target])
    edge_vec = np.array([r_edge_val, z_edge_val])

    diff_target = target_vec - axis_vec
    diff_edge = edge_vec - axis_vec
    edge_norm2 = diff_edge.dot(diff_edge)
    if edge_norm2 <= 0.0:
        return 0.5, theta_guess

    ratio = diff_target.dot(diff_edge) / edge_norm2
    ratio = max(0.0, min(1.0, ratio))
    s_guess = ratio ** 2
    return s_guess, theta_guess


def _solve_flux_coordinates(
    r_target: float,
    z_target: float,
    phi: float,
    s0: float,
    theta0: float,
    max_iter: int = 20,
    tol: float = 1.0e-6,
) -> Tuple[float, float] | Tuple[None, None]:
    s = float(min(max(s0, 0.0), 1.0))
    theta = theta0 % (2.0 * math.pi)

    for _ in range(max_iter):
        res = _VMEC_WRAPPERS.splint_vmec_data_wrapper(s, theta, phi)
        (
            _A_phi,
            _A_theta,
            _dA_phi_ds,
            _dA_theta_ds,
            _aiota,
            R_cm,
            Z_cm,
            _alam,
            dR_ds_cm,
            dR_dt_cm,
            _dR_dp_cm,
            dZ_ds_cm,
            dZ_dt_cm,
            _dZ_dp_cm,
            _dl_ds,
            _dl_dt,
            _dl_dp,
        ) = res

        R = float(R_cm) * CM_TO_M
        Z = float(Z_cm) * CM_TO_M
        F = np.array([R - r_target, Z - z_target])
        norm_f = np.linalg.norm(F)
        if norm_f < tol:
            return s, theta

        J = np.array(
            [
                [float(dR_ds_cm) * CM_TO_M, float(dR_dt_cm) * CM_TO_M],
                [float(dZ_ds_cm) * CM_TO_M, float(dZ_dt_cm) * CM_TO_M],
            ]
        )

        try:
            delta = np.linalg.solve(J, F)
        except np.linalg.LinAlgError:
            return (None, None)

        s -= delta[0]
        theta -= delta[1]
        s = min(max(s, -0.1), 1.1)
        theta %= 2.0 * math.pi

    return (None, None)


def _evaluate_bfield(
    s: float,
    theta: float,
    phi: float,
) -> Tuple[float, float, float]:
    BR, Bphi, BZ, _Bmag = _VMEC_WRAPPERS.vmec_field_cylindrical_wrapper(s, theta, phi)
    return float(BR), float(Bphi), float(BZ)

def field_from_vmec(
    wout_path: str | Path,
    nr: int = 64,
    nz: int = 64,
    nphi: int = 16,
    r_margin: float = 0.05,
    z_margin: float = 0.05,
    max_iter: int = 20,
    tol: float = 1.0e-6,
) -> B3DSField:
    wout = Path(wout_path)
    nfp, axis_r, axis_z, iota_profile = _load_vmec_metadata(wout)
    geom = VMECGeometry.from_file(str(wout))

    phi_samples, r_edge, z_edge = _prepare_outer_surface(geom, nphi)
    rmin = float(np.min(r_edge) * CM_TO_M) - r_margin
    rmax = float(np.max(r_edge) * CM_TO_M) + r_margin
    zmin = float(np.min(z_edge) * CM_TO_M) - z_margin
    zmax = float(np.max(z_edge) * CM_TO_M) + z_margin

    r_grid = np.linspace(rmin, rmax, nr)
    z_grid = np.linspace(zmin, zmax, nz)
    phi_grid_deg = np.linspace(0.0, 360.0 / nfp, nphi, endpoint=False)
    phi_grid_rad = np.deg2rad(phi_grid_deg)

    _magfie.init_vmec(str(wout), 5)

    br = np.zeros((nr, nphi, nz), dtype=float)
    bphi = np.zeros_like(br)
    bz = np.zeros_like(br)
    s_map = np.full((nr, nz), np.nan, dtype=float)

    prev_s = 0.5
    prev_theta = 0.0

    for iphi, phi in enumerate(phi_grid_rad):
        prev_s_row = 0.5
        prev_theta_row = 0.0
        for ir, r_val in enumerate(r_grid):
            prev_s_col = prev_s_row
            prev_theta_col = prev_theta_row
            for iz, z_val in enumerate(z_grid):
                s_guess, theta_guess = _initial_guess(
                    axis_r,
                    axis_z,
                    iphi,
                    phi_samples,
                    r_edge,
                    z_edge,
                    r_val,
                    z_val,
                )
                seed_s = prev_s_col if not math.isnan(prev_s_col) else s_guess
                seed_theta = prev_theta_col if not math.isnan(prev_theta_col) else theta_guess

                sol_s, sol_theta = _solve_flux_coordinates(
                    r_val,
                    z_val,
                    phi,
                    seed_s,
                    seed_theta,
                    max_iter=max_iter,
                    tol=tol,
                )

                if sol_s is None:
                    br[ir, iphi, iz] = 0.0
                    bphi[ir, iphi, iz] = 0.0
                    bz[ir, iphi, iz] = 0.0
                    prev_s_col = math.nan
                    prev_theta_col = math.nan
                    continue

                BR, BPHI, BZ = _evaluate_bfield(sol_s, sol_theta, phi)

                br[ir, iphi, iz] = BR
                bphi[ir, iphi, iz] = BPHI
                bz[ir, iphi, iz] = BZ
                prev_s_col = sol_s
                prev_theta_col = sol_theta

                if iphi == 0:
                    s_map[ir, iz] = sol_s

            prev_s_row = prev_s_col if not math.isnan(prev_s_col) else prev_s_row
            prev_theta_row = prev_theta_col if not math.isnan(prev_theta_col) else prev_theta_row

        prev_s = prev_s_row
        prev_theta = prev_theta_row

    psi_grid = np.zeros((nr, nz), dtype=float)
    stor_grid = np.linspace(0.0, 1.0, iota_profile.size)
    iota_abs = np.abs(iota_profile)
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (iota_abs[1:] + iota_abs[:-1]) * np.diff(stor_grid))])
    if integ[-1] > 0.0:
        spol_grid = integ / integ[-1]
    else:
        spol_grid = integ
    mask = np.isfinite(s_map)
    psi_grid[mask] = np.interp(s_map[mask], stor_grid, spol_grid)
    psi_grid[~mask] = 1.0

    return B3DSField(
        r_grid=r_grid,
        z_grid=z_grid,
        phi_grid_deg=phi_grid_deg,
        axis_r=axis_r,
        axis_z=axis_z,
        psi=psi_grid,
        psi0=0.0,
        psi1=1.0,
        br=br,
        bphi=bphi,
        bz=bz,
    )


def field_from_mgrid(
    mgrid_path: str | Path,
    wout_path: str | Path | None = None,
) -> B3DSField:
    mgrid = Path(mgrid_path)
    with netCDF4.Dataset(mgrid) as ds:
        br = np.array(ds["br_001"][...], dtype=float)
        bphi = np.array(ds["bp_001"][...], dtype=float)
        bz = np.array(ds["bz_001"][...], dtype=float)

        rmin = float(ds["rmin"][()])
        rmax = float(ds["rmax"][()])
        zmin = float(ds["zmin"][()])
        zmax = float(ds["zmax"][()])

        nr = br.shape[2]
        nz = br.shape[1]
        nphi = br.shape[0]

        r_grid = np.linspace(rmin, rmax, nr)
        z_grid = np.linspace(zmin, zmax, nz)

        if "phi" in ds.variables:
            phi_vals = np.array(ds["phi"][...], dtype=float)
            phi_grid_deg = np.rad2deg(phi_vals)
        else:
            nfp = int(ds.variables.get("nfp", np.array([1]))[0])
            phi_grid_deg = np.linspace(0.0, 360.0 / nfp, nphi, endpoint=False)

    br = np.transpose(br, (2, 0, 1))
    bphi = np.transpose(bphi, (2, 0, 1))
    bz = np.transpose(bz, (2, 0, 1))

    axis_r = 0.5 * (rmin + rmax)
    axis_z = 0.5 * (zmin + zmax)

    if wout_path is not None:
        _, _, _, iota_profile = _load_vmec_metadata(Path(wout_path))
        stor_grid = np.linspace(0.0, 1.0, iota_profile.size)
        iota_abs = np.abs(iota_profile)
        integ = np.concatenate([[0.0], np.cumsum(0.5 * (iota_abs[1:] + iota_abs[:-1]) * np.diff(stor_grid))])
        spol_grid = integ / integ[-1] if integ[-1] > 0 else integ
        s_guess = np.linspace(0.0, 1.0, nr)
        psi_line = np.interp(s_guess, stor_grid, spol_grid)
        psi = np.tile(psi_line[:, None], (1, nz))
    else:
        psi = np.zeros((nr, nz), dtype=float)

    return B3DSField(
        r_grid=r_grid,
        z_grid=z_grid,
        phi_grid_deg=np.asarray(phi_grid_deg, dtype=float),
        axis_r=axis_r,
        axis_z=axis_z,
        psi=psi,
        psi0=0.0,
        psi1=1.0,
        br=br,
        bphi=bphi,
        bz=bz,
    )


def write_b3ds_hdf5(path: str | Path, field: B3DSField, desc: str | None = None, activate: bool = True) -> str:
    """Write a :class:`B3DSField` to an ASCOT5-compatible HDF5 file."""

    path = Path(path)
    data = field.as_dict()
    path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(path, "a") as h5:
        grp = h5.require_group("bfield")
        existing = [name for name in grp.keys() if name.startswith("B_3DS")]
        postfix = int(time.time())
        while f"B_3DS_{postfix}" in existing:
            postfix += 1
        gname = f"B_3DS_{postfix}"
        g = grp.create_group(gname)

        g.create_dataset("b_rmin", data=(data["b_rmin"],), dtype="f8")
        g.create_dataset("b_rmax", data=(data["b_rmax"],), dtype="f8")
        g.create_dataset("b_nr", data=(data["b_nr"],), dtype="i4")
        g.create_dataset("b_zmin", data=(data["b_zmin"],), dtype="f8")
        g.create_dataset("b_zmax", data=(data["b_zmax"],), dtype="f8")
        g.create_dataset("b_nz", data=(data["b_nz"],), dtype="i4")
        g.create_dataset("b_phimin", data=(data["b_phimin"],), dtype="f8")
        g.create_dataset("b_phimax", data=(data["b_phimax"],), dtype="f8")
        g.create_dataset("b_nphi", data=(data["b_nphi"],), dtype="i4")
        g.create_dataset("axisr", data=(data["axisr"],), dtype="f8")
        g.create_dataset("axisz", data=(data["axisz"],), dtype="f8")
        g.create_dataset("psi0", data=(data["psi0"],), dtype="f8")
        g.create_dataset("psi1", data=(data["psi1"],), dtype="f8")
        g.create_dataset("psi", data=data["psi"], dtype="f8")
        br_si = np.transpose(field.br, (2, 1, 0)) * GAUSS_TO_TESLA
        bphi_si = np.transpose(field.bphi, (2, 1, 0)) * GAUSS_TO_TESLA
        bz_si = np.transpose(field.bz, (2, 1, 0)) * GAUSS_TO_TESLA

        g.create_dataset("br", data=br_si, dtype="f8")
        g.create_dataset("bphi", data=bphi_si, dtype="f8")
        g.create_dataset("bz", data=bz_si, dtype="f8")

        if desc is not None:
            g.attrs["desc"] = np.bytes_(desc)

        if activate:
            grp.attrs["active"] = np.bytes_(gname.split("_")[-1])

    return gname
