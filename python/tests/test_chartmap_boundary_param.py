from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


@pytest.mark.network
def test_chartmap_boundary_param_theta_vs_arc(tmp_path: Path) -> None:
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary
    from libneo.vmec import VMECGeometry

    wout_url = "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    theta_chart = tmp_path / "chartmap_theta.nc"
    arc_chart = tmp_path / "chartmap_arc.nc"

    write_chartmap_from_vmec_boundary(
        wout,
        theta_chart,
        nrho=33,
        ntheta=65,
        nzeta=33,
        s_boundary=1.0,
        boundary_offset=0.0,
        boundary_param="theta",
    )
    write_chartmap_from_vmec_boundary(
        wout,
        arc_chart,
        nrho=33,
        ntheta=65,
        nzeta=33,
        s_boundary=1.0,
        boundary_offset=0.0,
        boundary_param="arc",
    )

    geom = VMECGeometry.from_file(str(wout))

    def _chartmap_boundary(chartmap_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        import netCDF4

        with netCDF4.Dataset(chartmap_path, "r") as ds:
            theta = np.array(ds.variables["theta"][:], dtype=float)
            zeta = np.array(ds.variables["zeta"][:], dtype=float)
            x = np.array(ds.variables["x"][:], dtype=float)
            y = np.array(ds.variables["y"][:], dtype=float)
            z = np.array(ds.variables["z"][:], dtype=float)
            units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
            scale = 0.01 if units == "cm" else 1.0

        r = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
        z = z * scale
        ir = r.shape[2] - 1
        return theta, zeta, r[:, :, ir], z[:, :, ir]

    theta_vals, zeta_vals, R_theta, Z_theta = _chartmap_boundary(theta_chart)
    _, _, R_arc, Z_arc = _chartmap_boundary(arc_chart)

    zeta0 = float(zeta_vals[0])
    R_vmec, Z_vmec, _ = geom.boundary_rz(1.0, theta_vals, zeta0, boundary_offset=0.0, use_asym=True)

    err_theta = max(
        float(np.max(np.abs(R_theta[0, :] - R_vmec))),
        float(np.max(np.abs(Z_theta[0, :] - Z_vmec))),
    )
    err_arc = max(
        float(np.max(np.abs(R_arc[0, :] - R_vmec))),
        float(np.max(np.abs(Z_arc[0, :] - Z_vmec))),
    )

    assert err_theta < 2.0e-2
    assert err_arc > 5.0e-2
