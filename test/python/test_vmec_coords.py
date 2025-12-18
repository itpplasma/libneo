import os
import tempfile
from urllib.request import urlopen

import numpy as np
import pytest

from libneo.vmec import VMECGeometry, vmec_to_cylindrical
from libneo.vmec_to_efit import cfunct as ref_cfunct, sfunct as ref_sfunct


STELLOPT_WOUT_URL = (
    "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
)


def _download_file(url: str, dst: str) -> None:
    with urlopen(url, timeout=30) as resp, open(dst, "wb") as f:
        f.write(resp.read())


@pytest.mark.network
def test_vmec_axis_matches_wout_axis_coeffs():
    """
    Validate VMEC (R,Z) reconstruction against the explicit magnetic axis
    coefficients stored in the wout file.

    VMEC stores `zaxis_cs`, but the physical axis uses the VMEC phase convention.
    To keep the sign convention explicit and consistent with `sin(m*theta - n*zeta)`,
    write the axis as:

      Z_axis(phi) = sum_n zaxis_cs(n) * sin(-n*nfp*phi)

    (using sin(-x) = -sin(x)).
    """
    from netCDF4 import Dataset

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "wout.nc")
        _download_file(STELLOPT_WOUT_URL, path)

        geom = VMECGeometry.from_file(path)

        with Dataset(path) as nc:
            raxis_cc = np.array(nc.variables["raxis_cc"][:], dtype=float)
            zaxis_cs = np.array(nc.variables["zaxis_cs"][:], dtype=float)
            nfp = int(np.array(nc.variables["nfp"][:]))

        n = np.arange(raxis_cc.size, dtype=float)
        phis = np.array([0.2, 0.7, 1.1], dtype=float)
        for phi in phis:
            R_axis = np.sum(raxis_cc * np.cos(n * nfp * phi))
            Z_axis = np.sum(zaxis_cs * np.sin(-n * nfp * phi))

            R, Z, phi_out = geom.coords(s_index=0, theta=np.array([0.0]), zeta=float(phi), use_asym=True)
            assert phi_out == pytest.approx(phi)
            assert R.shape == Z.shape == (1,)
            assert np.allclose(R[0], R_axis, rtol=0.0, atol=1e-10)
            assert np.allclose(Z[0], Z_axis, rtol=0.0, atol=1e-10)


@pytest.mark.network
def test_vmec_to_cylindrical_matches_reference():
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "wout.nc")
        _download_file(STELLOPT_WOUT_URL, path)

        geom = VMECGeometry.from_file(path)
        ns = geom.rmnc.shape[1]
        s_idx = ns - 1  # outermost surface

        theta = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
        for zeta in (0.0, np.pi / 3.0):
            R, Z, phi = vmec_to_cylindrical(path, s_idx, theta, zeta, use_asym=True)

            # Reference via direct evaluation using vmec_to_efit helpers
            R_ref = ref_cfunct(theta, zeta, geom.rmnc, geom.xm, geom.xn)[s_idx, :]
            Z_ref = ref_sfunct(theta, zeta, geom.zmns, geom.xm, geom.xn)[s_idx, :]
            if geom.rmns is not None and geom.zmnc is not None:
                R_ref = R_ref + ref_sfunct(theta, zeta, geom.rmns, geom.xm, geom.xn)[s_idx, :]
                Z_ref = Z_ref + ref_cfunct(theta, zeta, geom.zmnc, geom.xm, geom.xn)[s_idx, :]

            assert np.allclose(R, R_ref, rtol=0.0, atol=1e-10)
            assert np.allclose(Z, Z_ref, rtol=0.0, atol=1e-10)
            assert phi == pytest.approx(zeta)


@pytest.mark.network
def test_vmec_coords_s_matches_linear_phi_interpolation():
    from netCDF4 import Dataset

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "wout.nc")
        _download_file(STELLOPT_WOUT_URL, path)

        geom = VMECGeometry.from_file(path)
        ns = int(geom.rmnc.shape[1])
        assert ns >= 2

        with Dataset(path) as nc:
            phi = np.array(nc.variables["phi"][:], dtype=float)
        phi_edge = float(phi[-1])

        theta = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
        zeta = 0.31

        s = 0.43
        phi_target = float(s) * phi_edge
        i1 = int(np.searchsorted(phi, phi_target, side="right"))
        i1 = max(1, min(ns - 1, i1))
        i0 = i1 - 1
        denom = float(phi[i1] - phi[i0])
        alpha = 0.0 if denom == 0.0 else float((phi_target - phi[i0]) / denom)

        R0, Z0, _ = geom.coords(i0, theta, zeta, use_asym=True)
        R1, Z1, _ = geom.coords(i1, theta, zeta, use_asym=True)
        R_ref = (1.0 - alpha) * R0 + alpha * R1
        Z_ref = (1.0 - alpha) * Z0 + alpha * Z1

        R_s, Z_s, zeta_out = geom.coords_s(s, theta, zeta, use_asym=True)
        assert zeta_out == pytest.approx(zeta)
        assert np.allclose(R_s, R_ref, rtol=0.0, atol=1e-10)
        assert np.allclose(Z_s, Z_ref, rtol=0.0, atol=1e-10)


@pytest.mark.network
def test_vmec_plot_surfaces_visual_check(tmp_path):
    """
    Visual check: plot several nested flux surface cross-sections at
    multiple zeta values in the RZ plane and save to file.
    """
    import matplotlib.pyplot as plt
    from pathlib import Path

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "wout.nc")
        _download_file(STELLOPT_WOUT_URL, path)

        geom = VMECGeometry.from_file(path)
        ns = geom.rmnc.shape[1]
        # Pick ~5 surfaces spanning core to edge
        s_indices = np.linspace(0, ns - 1, 5, dtype=int)
        thetas = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
        zetas = [0.0, np.pi / 4.0, np.pi / 2.0]

        fig, axes = plt.subplots(1, len(zetas), figsize=(12, 4), constrained_layout=True)
        if len(zetas) == 1:
            axes = [axes]
        spans_RZ = []
        for ax, zeta in zip(axes, zetas):
            for s_idx in s_indices:
                R, Z, _ = vmec_to_cylindrical(path, s_idx, thetas, zeta, use_asym=True)
                ax.plot(R, Z, lw=1)
                # Basic sanity: shapes and finiteness
                assert R.shape == Z.shape == thetas.shape
                assert np.isfinite(R).all() and np.isfinite(Z).all()
                # Track spans for monotonic growth check
                spans_RZ.append((zeta, s_idx, R.max() - R.min(), Z.max() - Z.min()))
            ax.set_title(f"zeta={zeta:.2f} rad")
            ax.set_xlabel("R [m]")
            ax.set_ylabel("Z [m]")
            ax.axis("equal")

        # Save to build directory for artifact collection
        build_artifact = Path("build/test/python/vmec_surfaces.png")
        build_artifact.parent.mkdir(parents=True, exist_ok=True)
        outfile = build_artifact
        fig.savefig(build_artifact)
        plt.close(fig)
        assert outfile.exists()
        assert outfile.stat().st_size > 0

        # Sanity: for a fixed zeta, outer surfaces should have larger span
        # Sort by s_idx for first zeta only
        z0 = zetas[0]
        spans_z0 = [(s, dR, dZ) for z, s, dR, dZ in spans_RZ if np.isclose(z, z0)]
        spans_z0.sort(key=lambda t: t[0])
        # Require strictly increasing R-span for at least consecutive pairs
        if len(spans_z0) >= 2:
            dRs = [dR for _, dR, _ in [(s, dR, dZ) for s, dR, dZ in spans_z0]]
            assert all(dRs[i] < dRs[i+1] for i in range(len(dRs)-1))

        # Sanity: periodicity at theta=0 and 2π (evaluate explicit closure)
        thetas_closed = np.concatenate([thetas, [2.0 * np.pi]])
        for s_idx in s_indices[:2]:  # a couple surfaces suffice
            Rc, Zc, _ = vmec_to_cylindrical(path, s_idx, thetas_closed, z0, use_asym=True)
            assert abs(Rc[0] - Rc[-1]) < 1e-8
            assert abs(Zc[0] - Zc[-1]) < 1e-8
