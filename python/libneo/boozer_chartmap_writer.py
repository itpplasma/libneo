"""Write an extended Boozer chartmap NetCDF file.

The extended chartmap is the on-disk format SIMPLE reads without any VMEC,
booz_xform or GVEC library at runtime. It stores the Boozer-coordinate
geometry (x, y, z on a rho/theta/zeta grid), the magnetic field strength
Bmod on the same grid, and the radial surface functions A_phi, B_theta,
B_phi. The format is defined by libneo; converter front ends (booz_xform,
GVEC) fill the arrays and call write_boozer_chartmap.

Units: inputs are SI (T, m, T*m, T*m^2); the chartmap stores CGS-Gaussian
(G, cm, G*cm, G*cm^2). The SI->CGS factors below are the single source of
truth shared by all front ends.

Handedness: SIMPLE uses left-handed Boozer coordinates. apply_left_handed_flip
negates the toroidal surface functions (B_phi, A_phi) for a toroidal-angle
flip, matching the convention in SIMPLE's converters.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# SI -> CGS-Gaussian unit conversions (single source of truth).
TESLA_TO_GAUSS = 1.0e4  # T -> G
METER_TO_CM = 1.0e2  # m -> cm
TESLA_METER_TO_GAUSS_CM = 1.0e6  # T*m -> G*cm
TESLA_METER2_TO_GAUSS_CM2 = 1.0e8  # T*m^2 -> G*cm^2


def apply_left_handed_flip(B_phi, A_phi):
    """Negate toroidal surface functions for a left-handed toroidal flip.

    Right-handed (VMEC/GVEC) -> left-handed (SIMPLE) coordinates. Returns the
    flipped (B_phi, A_phi). The poloidal-flip variant (negating B_theta and
    the toroidal flux) is left to the caller because it depends on which angle
    the geometry grid was built against.
    """
    return -np.asarray(B_phi), -np.asarray(A_phi)


def write_boozer_chartmap(
    path,
    *,
    rho,
    s,
    theta,
    zeta,
    x,
    y,
    z,
    A_phi,
    B_theta,
    B_phi,
    Bmod,
    num_field_periods,
    torflux,
    **attrs,
):
    """Write an extended Boozer chartmap NetCDF file.

    All field/geometry inputs are in CGS-Gaussian units already (callers use
    the module constants to convert from SI). The (rho, theta, zeta) arrays
    are 1D abscissae; x, y, z and Bmod are (n_rho, n_theta, n_zeta) and are
    transposed to NetCDF (zeta, theta, rho) on write so NF90 reads them back
    as (rho, theta, zeta).

    Parameters
    ----------
    path : path-like
        Output .nc file.
    rho, s : 1D arrays, length n_rho
        Radial abscissae (rho = effective radius, s = normalized tor. flux).
    theta, zeta : 1D arrays
        Poloidal and toroidal Boozer angles (zeta over one field period).
    x, y, z : (n_rho, n_theta, n_zeta) arrays
        Cartesian Boozer geometry in cm.
    A_phi : 1D array, length n_rho
        Toroidal vector potential surface function on s, in G*cm^2.
    B_theta, B_phi : 1D arrays, length n_rho
        Covariant B surface functions on rho, in G*cm.
    Bmod : (n_rho, n_theta, n_zeta) array
        Field strength in G.
    num_field_periods : int
        Number of field periods.
    torflux : float
        Toroidal flux (A_theta at the edge) in G*cm^2.
    **attrs
        Extra global attributes passed through verbatim (provenance, e.g.
        booz2chartmap_source or gvec2chartmap_flip).
    """
    from netCDF4 import Dataset

    rho = np.asarray(rho, dtype=float)
    s = np.asarray(s, dtype=float)
    theta = np.asarray(theta, dtype=float)
    zeta = np.asarray(zeta, dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    A_phi = np.asarray(A_phi, dtype=float)
    B_theta = np.asarray(B_theta, dtype=float)
    B_phi = np.asarray(B_phi, dtype=float)
    Bmod = np.asarray(Bmod, dtype=float)

    n_rho = rho.size
    n_theta = theta.size
    n_zeta = zeta.size
    shape = (n_rho, n_theta, n_zeta)
    for name, arr in [("x", x), ("y", y), ("z", z), ("Bmod", Bmod)]:
        if arr.shape != shape:
            raise ValueError(
                f"{name} has shape {arr.shape}, expected {shape} "
                "(n_rho, n_theta, n_zeta)"
            )
    if s.size != n_rho:
        raise ValueError(f"s has length {s.size}, expected {n_rho}")
    if A_phi.size != n_rho:
        raise ValueError(f"A_phi has length {A_phi.size}, expected {n_rho}")
    for name, arr in [("B_theta", B_theta), ("B_phi", B_phi)]:
        if arr.size != n_rho:
            raise ValueError(f"{name} has length {arr.size}, expected {n_rho}")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    ds = Dataset(path, "w", format="NETCDF4")
    try:
        ds.createDimension("rho", n_rho)
        ds.createDimension("s", n_rho)
        ds.createDimension("theta", n_theta)
        ds.createDimension("zeta", n_zeta)

        v = ds.createVariable("rho", "f8", ("rho",))
        v[:] = rho
        v = ds.createVariable("s", "f8", ("s",))
        v[:] = s
        v = ds.createVariable("theta", "f8", ("theta",))
        v[:] = theta
        v = ds.createVariable("zeta", "f8", ("zeta",))
        v[:] = zeta

        # NetCDF dims (zeta, theta, rho) so NF90 reads them as (rho, theta, zeta).
        for name, arr in [("x", x), ("y", y), ("z", z)]:
            v = ds.createVariable(name, "f8", ("zeta", "theta", "rho"))
            v[:] = np.transpose(arr, (2, 1, 0))
            v.units = "cm"

        v = ds.createVariable("A_phi", "f8", ("s",))
        v[:] = A_phi
        v.radial_abscissa = "s"
        v = ds.createVariable("B_theta", "f8", ("rho",))
        v[:] = B_theta
        v = ds.createVariable("B_phi", "f8", ("rho",))
        v[:] = B_phi
        v = ds.createVariable("Bmod", "f8", ("zeta", "theta", "rho"))
        v[:] = np.transpose(Bmod, (2, 1, 0))

        v = ds.createVariable("num_field_periods", "i4")
        v.assignValue(np.int32(num_field_periods))

        ds.rho_convention = attrs.pop("rho_convention", "rho_tor")
        ds.zeta_convention = attrs.pop("zeta_convention", "boozer")
        ds.rho_lcfs = float(attrs.pop("rho_lcfs", float(rho[-1])))
        ds.boozer_field = np.int32(attrs.pop("boozer_field", 1))
        ds.torflux = float(torflux)

        for key, value in attrs.items():
            ds.setncattr(key, value)
    finally:
        ds.close()
