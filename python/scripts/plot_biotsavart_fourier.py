#!/usr/bin/env python3
"""Quick interactive plot comparing Fourier-mode magnitudes from two coil_tools outputs."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from libneo.biotsavart_fourier import (
    read_Bnvac_fourier,
    read_Anvac_fourier,
    gauged_Anvac_from_Bnvac,
    gauge_Anvac,
    spline_gauged_Anvac,
    field_divfree,
)


def main() -> None:
    ref_path = Path(input("Path to coil_tools Fourier HDF5 file: ").strip())
    test_path = Path(input("Path to vector-potential NetCDF file: ").strip())

    ref_grid, BnR_ref, Bnphi_ref, BnZ_ref = read_Bnvac_fourier(str(ref_path))
    test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
        str(test_path)
    )

    ref_spl = spline_gauged_Anvac(
        ref_grid, *gauged_Anvac_from_Bnvac(ref_grid, BnR_ref, Bnphi_ref, BnZ_ref)
    )
    test_spl = spline_gauged_Anvac(
        test_grid, *gauge_Anvac(test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)
    )

    if not np.allclose(ref_grid.R, test_grid.R) or not np.allclose(
        ref_grid.Z, test_grid.Z
    ):
        raise ValueError(
            "Reference and test grids differ; please resample before plotting"
        )

    R = np.linspace(ref_grid.R_min, ref_grid.R_max, 2 * ref_grid.nR - 1)
    Z = np.linspace(ref_grid.Z_min, ref_grid.Z_max, 2 * ref_grid.nZ - 1)

    ncoil = len(ref_spl["AnR_Re"])
    log_bn2 = np.empty((2, ncoil, R.size, Z.size), dtype=float)

    for idx, spl in enumerate((ref_spl, test_spl)):
        BnR, Bnphi, BnZ = field_divfree(spl, R, Z)
        magnitude_sq = (
            BnR * np.conj(BnR) + Bnphi * np.conj(Bnphi) + BnZ * np.conj(BnZ)
        ).real
        log_bn2[idx] = np.log10(np.maximum(magnitude_sq, 1e-300))

    norm = Normalize(vmin=np.min(log_bn2), vmax=np.max(log_bn2))
    fig, axs = plt.subplots(ncoil, 2, figsize=(6, 3 * ncoil), layout="constrained")

    for coil in range(ncoil):
        for column, label in enumerate(("Fourier", "vector")):
            ax = axs[coil, column] if ncoil > 1 else axs[column]
            im = ax.imshow(
                log_bn2[column][coil].T,
                origin="lower",
                cmap="magma",
                extent=[R[0], R[-1], Z[0], Z[-1]],
                norm=norm,
                aspect="auto",
                interpolation="bilinear",
            )
            ax.set_title(f"Coil {coil + 1}: {label}")
            ax.set_xlabel("R [cm]")
            ax.set_ylabel("Z [cm]")

    cbar = fig.colorbar(im, ax=axs, location="bottom", fraction=0.05, pad=0.08)
    cbar.set_label(r"$\log_{10} |\vec{B}_n|^2$")
    plt.show()


if __name__ == "__main__":
    main()
