class grid_t:
    def __init__(self):
        from numpy import empty

        self.nR = 0
        self.nZ = 0
        self.R_min = 0.0
        self.R_max = 0.0
        self.Z_min = 0.0
        self.Z_max = 0.0
        self.R = empty((0,))
        self.Z = empty((0,))


def read_Bnvac_fourier(field_file='AUG_B_coils.h5', ntor=2):
    from os import path, strerror
    from errno import ENOENT
    from numpy import empty, linspace
    import h5py

    if (not path.exists(field_file)):
        raise FileNotFoundError(ENOENT, strerror(ENOENT), field_file)
    h5py.get_config().complex_names = ('real', 'imag')
    grid = grid_t()
    with h5py.File(field_file, 'r') as file:
        ncoil = file[f"/ntor_{ntor:02}/ncoil"][()]
        grid.nR = file[f"/ntor_{ntor:02}/nR"][()]
        grid.nZ = file[f"/ntor_{ntor:02}/nZ"][()]
        grid.R_min = file[f"/ntor_{ntor:02}/R_min"][()]
        grid.R_max = file[f"/ntor_{ntor:02}/R_max"][()]
        grid.Z_min = file[f"/ntor_{ntor:02}/Z_min"][()]
        grid.Z_max = file[f"/ntor_{ntor:02}/Z_max"][()]
        grid.R = linspace(grid.R_min, grid.R_max, grid.nR)
        grid.Z = linspace(grid.Z_min, grid.Z_max, grid.nZ)
        BnR = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
        Bnphi = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
        BnZ = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
        for kcoil in range(1, ncoil + 1):
            BnR[kcoil - 1, :, :] = file[f"/ntor_{ntor:02}/coil_{kcoil:02}/Bn_R"][()].T
            Bnphi[kcoil - 1, :, :] = file[f"/ntor_{ntor:02}/coil_{kcoil:02}/Bn_phi"][()].T
            BnZ[kcoil - 1, :, :] = file[f"/ntor_{ntor:02}/coil_{kcoil:02}/Bn_Z"][()].T
    return grid, BnR, Bnphi, BnZ


def read_Anvac_fourier(field_file='AUG_B_coils.nc', ntor=2):
    from os import path, strerror
    from errno import ENOENT
    from numpy import array, empty
    import netCDF4

    if (not path.exists(field_file)):
        raise FileNotFoundError(ENOENT, strerror(ENOENT), field_file)
    grid = grid_t()
    rootgrp = netCDF4.Dataset(field_file, 'r')
    ncoil = len(rootgrp.dimensions['coil_number'])
    grid.R = array(rootgrp['R'])
    grid.Z = array(rootgrp['Z'])
    grid.nR = grid.R.size
    grid.nZ = grid.Z.size
    grid.R_min = grid.R[0]
    grid.R_max = grid.R[-1]
    grid.Z_min = grid.Z[0]
    grid.Z_max = grid.Z[-1]
    AnR = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
    Anphi = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
    AnZ = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
    dAnphi_dR = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
    dAnphi_dZ = empty((ncoil, grid.nR, grid.nZ), dtype=complex)
    for kcoil in range(ncoil):
        AnR[kcoil, :, :].real = rootgrp['AnR_real'][kcoil, :, :, ntor].T
        AnR[kcoil, :, :].imag = rootgrp['AnR_imag'][kcoil, :, :, ntor].T
        Anphi[kcoil, :, :].real = rootgrp['Anphi_real'][kcoil, :, :, ntor].T
        Anphi[kcoil, :, :].imag = rootgrp['Anphi_imag'][kcoil, :, :, ntor].T
        AnZ[kcoil, :, :].real = rootgrp['AnZ_real'][kcoil, :, :, ntor].T
        AnZ[kcoil, :, :].imag = rootgrp['AnZ_imag'][kcoil, :, :, ntor].T
        dAnphi_dR[kcoil, :, :].real = rootgrp['dAnphi_dR_real'][kcoil, :, :, ntor].T
        dAnphi_dR[kcoil, :, :].imag = rootgrp['dAnphi_dR_imag'][kcoil, :, :, ntor].T
        dAnphi_dZ[kcoil, :, :].real = rootgrp['dAnphi_dZ_real'][kcoil, :, :, ntor].T
        dAnphi_dZ[kcoil, :, :].imag = rootgrp['dAnphi_dZ_imag'][kcoil, :, :, ntor].T
    rootgrp.close()
    return grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ


def gauged_Anvac_from_Bnvac(grid, BnR, BnZ, ntor=2):
    from numpy import newaxis

    gauged_AnR = 1j / ntor * grid.R[newaxis, :, newaxis] * BnZ
    gauged_AnZ = -1j / ntor * grid.R[newaxis, :, newaxis] * BnR
    return gauged_AnR, gauged_AnZ


def gauge_Anvac(grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=2):
    from numpy import newaxis

    if ntor == 0:
        # For axisymmetric case, no gauge transformation needed
        return AnR, AnZ

    gauged_AnR = AnR + 1j / ntor * grid.R[newaxis, :, newaxis] * dAnphi_dR + 1j / ntor * Anphi
    gauged_AnZ = AnZ + 1j / ntor * grid.R[newaxis, :, newaxis] * dAnphi_dZ
    return gauged_AnR, gauged_AnZ


def spline_gauged_Anvac(grid, gauged_AnR, gauged_AnZ, ntor=2, Anphi=None):
    from scipy.interpolate import RectBivariateSpline

    ncoil = gauged_AnR.shape[0]
    spl = {
        'AnR_Re': [None] * ncoil,
        'AnR_Im': [None] * ncoil,
        'AnZ_Re': [None] * ncoil,
        'AnZ_Im': [None] * ncoil,
        'ntor': ntor,
    }
    for kcoil in range(ncoil):
        spl['AnR_Re'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnR[kcoil].real, kx=5, ky=5)
        spl['AnR_Im'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnR[kcoil].imag, kx=5, ky=5)
        spl['AnZ_Re'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnZ[kcoil].real, kx=5, ky=5)
        spl['AnZ_Im'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnZ[kcoil].imag, kx=5, ky=5)

    # For ntor=0, also spline Aphi since we need it for B = curl(A)
    if ntor == 0:
        if Anphi is None:
            raise ValueError("For ntor=0, Anphi must be provided to spline_gauged_Anvac")
        spl['Anphi_Re'] = [None] * ncoil
        spl['Anphi_Im'] = [None] * ncoil
        for kcoil in range(ncoil):
            spl['Anphi_Re'][kcoil] = RectBivariateSpline(grid.R, grid.Z, Anphi[kcoil].real, kx=5, ky=5)
            spl['Anphi_Im'][kcoil] = RectBivariateSpline(grid.R, grid.Z, Anphi[kcoil].imag, kx=5, ky=5)

    return spl


def field_divfree(spl, R, Z, ntor=None):
    """Evaluate vector potential splines and return Fourier amplitude of magnetic field.

    Keyword arguments:
    spl -- dict of splines from spline_gauged_Anvac
    R -- array of radii in cylindrical coordinates
    Z -- array of altitudes in cylindrical coordinates
    ntor -- toroidal mode number (if None, read from spl dict)

    Return values:
    BnR -- R component of the magnetic field's Fourier mode
    Bnphi -- phi component of the magnetic field's Fourier mode
    BnZ -- Z component of the magnetic field's Fourier mode

    The returned arrays have shape (ncoil, nR, nZ), where nR and nZ
    are the number of elements in the arrays R and Z, respectively.
    ncoil is the number of coils taken from the input splines.
    """

    from numpy import atleast_1d, empty, newaxis, squeeze

    # Auto-detect ntor from spline dict if not provided
    if ntor is None:
        if 'ntor' in spl:
            ntor = spl['ntor']
        else:
            ntor = 2  # Default for backward compatibility

    R = atleast_1d(R).ravel()
    Z = atleast_1d(Z).ravel()
    nR = R.size
    nZ = Z.size
    ncoil = len(spl['AnR_Re'])
    AnR = empty((ncoil, nR, nZ), dtype=complex)
    AnZ = empty((ncoil, nR, nZ), dtype=complex)
    dAnR_dZ = empty((ncoil, nR, nZ), dtype=complex)
    dAnZ_dR = empty((ncoil, nR, nZ), dtype=complex)
    for kcoil in range(ncoil):
        AnR[kcoil, :, :] = spl['AnR_Re'][kcoil](R, Z) + 1j * spl['AnR_Im'][kcoil](R, Z)
        AnZ[kcoil, :, :] = spl['AnZ_Re'][kcoil](R, Z) + 1j * spl['AnZ_Im'][kcoil](R, Z)
        dAnR_dZ[kcoil, :, :] = spl['AnR_Re'][kcoil](R, Z, dy=1) + 1j * spl['AnR_Im'][kcoil](R, Z, dy=1)
        dAnZ_dR[kcoil, :, :] = spl['AnZ_Re'][kcoil](R, Z, dx=1) + 1j * spl['AnZ_Im'][kcoil](R, Z, dx=1)

    if ntor == 0:
        # For axisymmetric case (ntor=0), use curl in cylindrical coords:
        # B_R = -∂A_φ/∂Z
        # B_φ = ∂A_R/∂Z - ∂A_Z/∂R
        # B_Z = ∂A_φ/∂R + A_φ/R
        if 'Anphi_Re' not in spl:
            raise ValueError("For ntor=0, Aphi splines must be provided")

        Anphi = empty((ncoil, nR, nZ), dtype=complex)
        dAnphi_dR = empty((ncoil, nR, nZ), dtype=complex)
        dAnphi_dZ = empty((ncoil, nR, nZ), dtype=complex)
        for kcoil in range(ncoil):
            Anphi[kcoil, :, :] = spl['Anphi_Re'][kcoil](R, Z) + 1j * spl['Anphi_Im'][kcoil](R, Z)
            dAnphi_dR[kcoil, :, :] = spl['Anphi_Re'][kcoil](R, Z, dx=1) + 1j * spl['Anphi_Im'][kcoil](R, Z, dx=1)
            dAnphi_dZ[kcoil, :, :] = spl['Anphi_Re'][kcoil](R, Z, dy=1) + 1j * spl['Anphi_Im'][kcoil](R, Z, dy=1)

        BnR = -dAnphi_dZ
        Bnphi = dAnR_dZ - dAnZ_dR
        BnZ = dAnphi_dR + Anphi / R[newaxis, :, newaxis]
    else:
        # Original formulas for n ≠ 0
        BnR = 1j * ntor * AnZ / R[newaxis, :, newaxis]
        Bnphi = dAnR_dZ - dAnZ_dR
        BnZ = -1j * ntor * AnR / R[newaxis, :, newaxis]

    return squeeze(BnR), squeeze(Bnphi), squeeze(BnZ)
