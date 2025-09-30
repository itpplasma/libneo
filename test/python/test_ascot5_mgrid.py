import numpy as np
import matplotlib
matplotlib.use("Agg")  # pragma: no cover
import matplotlib.pyplot as plt
import netCDF4

from libneo.ascot5 import field_from_mgrid, write_b3ds_hdf5


def _create_synthetic_mgrid(path):
    nr, nz, nphi = 6, 5, 4
    rmin, rmax = 4.0, 5.5
    zmin, zmax = -0.5, 0.5

    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("phi", nphi)
        ds.createDimension("z", nz)
        ds.createDimension("r", nr)

        var_br = ds.createVariable("br_001", "f8", ("phi", "z", "r"))
        var_bphi = ds.createVariable("bp_001", "f8", ("phi", "z", "r"))
        var_bz = ds.createVariable("bz_001", "f8", ("phi", "z", "r"))

        phi_idx = np.arange(nphi)[:, None, None]
        z_idx = np.arange(nz)[None, :, None]
        r_idx = np.arange(nr)[None, None, :]

        var_br[:] = 0.002 * np.sin(2 * np.pi * phi_idx / nphi) * (z_idx - nz / 2)
        var_bphi[:] = 0.1 + 0.001 * (r_idx - nr / 2)
        var_bz[:] = 0.002 * np.cos(2 * np.pi * phi_idx / nphi) * (z_idx - nz / 2)

        ds.createVariable("rmin", "f8")[:] = rmin
        ds.createVariable("rmax", "f8")[:] = rmax
        ds.createVariable("zmin", "f8")[:] = zmin
        ds.createVariable("zmax", "f8")[:] = zmax
        ds.createVariable("nfp", "i4")[:] = 1

        phi_vals = ds.createVariable("phi", "f8", ("phi",))
        phi_vals[:] = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)


def test_field_from_mgrid(tmp_path):
    mgrid_path = tmp_path / "synthetic_mgrid.nc"
    _create_synthetic_mgrid(mgrid_path)

    field = field_from_mgrid(mgrid_path)
    assert field.br.shape == (6, 4, 5)
    assert np.isfinite(field.bphi).all()

    output = tmp_path / "bfield_mgrid.h5"
    gname = write_b3ds_hdf5(output, field, desc="mgrid-test")
    assert output.exists()
    assert gname.startswith("B_3DS")

    bmag = np.sqrt(field.br[:, 0, :] ** 2 + field.bphi[:, 0, :] ** 2 + field.bz[:, 0, :] ** 2).T
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.pcolormesh(field.r_grid, field.z_grid, bmag, shading="auto")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    fig.colorbar(im, ax=ax, label="|B| [T]")
    png_path = tmp_path / "mgrid_bfield.png"
    fig.savefig(png_path)
    plt.close(fig)
    assert png_path.exists()
