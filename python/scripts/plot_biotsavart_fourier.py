# %%
from numpy import amin, amax, arange, arctan2, conj, empty, linspace, log10, pi, sqrt
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from libneo.biotsavart_fourier import *

# %%
ref = input('Paste the path to the HDF5 file')
test = input('Paste the path to the NetCDF4 file')

# %%
ref_grid, BnR, Bnphi, BnZ = read_Bnvac_fourier(ref)
test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(test)

# %%
ref_spl = spline_gauged_Anvac(ref_grid,
    *gauged_Anvac_from_Bnvac(ref_grid, BnR, Bnphi, BnZ))
test_spl = spline_gauged_Anvac(test_grid,
    *gauge_Anvac(test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ))

# %%
assert(all(ref_grid.R == test_grid.R))
assert(all(ref_grid.Z == test_grid.Z))
# double the resolution
R = linspace(ref_grid.R_min, ref_grid.R_max, 2 * ref_grid.nR - 1)
Z = linspace(ref_grid.Z_min, ref_grid.Z_max, 2 * ref_grid.nZ - 1)


# %%
ncoil = len(ref_spl['AnR_Re'])
log_Bn2 = empty((2, ncoil, R.size, Z.size))
arg_Bn2 = empty((2, ncoil, R.size, Z.size))
for k, spl in enumerate((ref_spl, test_spl)):
    BnR, Bnphi, BnZ = field_divfree(spl, R, Z)
    log_Bn2[k] = log10(BnR * conj(BnR) + 
                       Bnphi * conj(Bnphi) + 
                       BnZ * conj(BnZ)).real
    arg_Bn2[k] = arctan2(
        sqrt(BnR.imag ** 2 + Bnphi.imag ** 2 + BnZ.imag ** 2),
        sqrt(BnR.real ** 2 + Bnphi.real ** 2 + BnZ.real ** 2)
    )

# %%
norm = Normalize(vmin=amin(log_Bn2), vmax=amax(log_Bn2))
fig = plt.figure(layout='constrained', figsize=(5, 3 * ncoil))
axs = fig.subplots(ncoil, 2)
for kcoil in range(ncoil):
    for k in range(2):
        im = axs[kcoil, k].imshow(log_Bn2[k][kcoil].T, origin='lower', cmap='magma',
                                  extent=[R[0], R[-1], Z[0], Z[-1]], norm=norm)
        axs[kcoil, k].set_title(f"coil {kcoil + 1}, {'reference' if k == 0 else 'test'}")
cbar = fig.colorbar(im, ax=axs, location='bottom')
cbar.set_label(r'$\log_{10} |\vec{B}_{n}|^{2}$')
plt.show()
