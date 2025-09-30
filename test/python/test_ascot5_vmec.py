import os
import shutil
import tempfile
from pathlib import Path
from urllib.request import urlopen

import matplotlib

matplotlib.use("Agg")  # pragma: no cover
import matplotlib.pyplot as plt
import numpy as np
import pytest

from libneo.ascot5 import GAUSS_TO_TESLA, field_from_vmec, write_b3ds_hdf5


STELLOPT_WOUT_URL = (
    "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
)


def _store_artifacts(*paths: Path) -> None:
    target_env = os.environ.get("PYTEST_ARTIFACTS")
    target = Path(target_env) if target_env else Path("pytest_artifacts")
    os.makedirs(target, exist_ok=True)
    for path in paths:
        shutil.copy(path, Path(target) / path.name)


@pytest.mark.network
def test_field_from_vmec_generates_b3ds(tmp_path):
    with tempfile.TemporaryDirectory() as td:
        wout_path = os.path.join(td, "wout.nc")
        with urlopen(STELLOPT_WOUT_URL, timeout=30) as resp, open(wout_path, "wb") as f:
            f.write(resp.read())

        field = field_from_vmec(
            wout_path,
            nr=24,
            nz=24,
            nphi=12,
            max_iter=30,
            tol=5.0e-6,
        )

    assert field.br.shape == (24, 12, 24)
    assert np.isfinite(field.br).any()
    bmag_full = np.sqrt(field.br**2 + field.bphi**2 + field.bz**2)
    finite_mask = np.isfinite(bmag_full)
    assert finite_mask.any(), "no finite magnetic field samples"
    assert finite_mask.mean() > 0.3
    nonzero_fraction = np.count_nonzero(bmag_full[finite_mask]) / finite_mask.sum()
    assert nonzero_fraction > 0.15

    output = tmp_path / "bfield_vmec.h5"
    gname = write_b3ds_hdf5(output, field, desc="vmec-test")
    assert output.exists()
    assert gname.startswith("B_3DS")

    # Produce diagnostic PNG
    bmag = (
        np.sqrt(
            field.br[:, 0, :] ** 2
            + field.bphi[:, 0, :] ** 2
            + field.bz[:, 0, :] ** 2
        ).T
        * GAUSS_TO_TESLA
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.pcolormesh(field.r_grid, field.z_grid, bmag, shading="auto")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    fig.colorbar(im, ax=ax, label="|B| [T]")
    png_path = tmp_path / "vmec_bfield.png"
    fig.savefig(png_path)
    plt.show()
    plt.close(fig)
    assert png_path.exists()
    _store_artifacts(output, png_path)
