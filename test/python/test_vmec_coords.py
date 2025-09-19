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
def test_vmec_plot_surfaces_visual_check(tmp_path):
    """
    Visual check: plot several nested flux surface cross-sections at
    multiple zeta values in the RZ plane and save to file.
    """
    import matplotlib.pyplot as plt

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
        outfile = tmp_path / "vmec_surfaces.png"
        fig.savefig(outfile)
        # Show interactively; in headless backends this is a no-op
        plt.show()
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
