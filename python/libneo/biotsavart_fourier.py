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


def read_Bnvac_fourier_all(field_file='AUG_B_coils.h5'):
    """Load all stored Fourier modes from a coil_tools Bnvac file.

    Returns
    -------
    grid : grid_t
        Cylindrical R/Z grid description.
    mode_numbers : ndarray
        1-D array of non-negative toroidal mode numbers present in the file.
    BnR, Bnphi, BnZ : ndarray
        Arrays with shape (nmodes, ncoil, nR, nZ) containing the complex
        Fourier amplitudes for each coil and mode.
    """

    from os import path, strerror
    from errno import ENOENT
    import h5py
    from numpy import array, linspace, empty, empty_like, asarray

    if not path.exists(field_file):
        raise FileNotFoundError(ENOENT, strerror(ENOENT), field_file)

    h5py.get_config().complex_names = ('real', 'imag')
    grid = grid_t()

    with h5py.File(field_file, 'r') as file:
        mode_numbers = sorted(
            int(key.split('_')[1])
            for key in file.keys()
            if key.startswith('ntor_')
        )
        if not mode_numbers:
            raise ValueError(f"No ntor groups found in {field_file}")

        first_group = file[f"/ntor_{mode_numbers[0]:02}"]
        ncoil = first_group['ncoil'][()]
        grid.nR = first_group['nR'][()]
        grid.nZ = first_group['nZ'][()]
        grid.R_min = first_group['R_min'][()]
        grid.R_max = first_group['R_max'][()]
        grid.Z_min = first_group['Z_min'][()]
        grid.Z_max = first_group['Z_max'][()]
        grid.R = linspace(grid.R_min, grid.R_max, grid.nR)
        grid.Z = linspace(grid.Z_min, grid.Z_max, grid.nZ)

        nmodes = len(mode_numbers)
        BnR = empty((nmodes, ncoil, grid.nR, grid.nZ), dtype=complex)
        Bnphi = empty_like(BnR)
        BnZ = empty_like(BnR)

        for idx, ntor in enumerate(mode_numbers):
            group_name = f"/ntor_{ntor:02}"
            g = file[group_name]
            if g['ncoil'][()] != ncoil:
                raise ValueError(f"Inconsistent coil count in group {group_name}")

            for kcoil in range(1, ncoil + 1):
                coil_group = f"{group_name}/coil_{kcoil:02}"
                BnR[idx, kcoil - 1] = file[coil_group + '/Bn_R'][()].T
                Bnphi[idx, kcoil - 1] = file[coil_group + '/Bn_phi'][()].T
                BnZ[idx, kcoil - 1] = file[coil_group + '/Bn_Z'][()].T

        # Convert structured dtype to complex if necessary
        if BnR.dtype.fields is not None:
            def _structured_to_complex(arr):
                return arr['real'] + 1j * arr['imag']

            BnR = _structured_to_complex(BnR)
            Bnphi = _structured_to_complex(Bnphi)
            BnZ = _structured_to_complex(BnZ)

    return grid, asarray(mode_numbers, dtype=int), BnR, Bnphi, BnZ


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


def read_Anvac_fourier_all(field_file='AUG_B_coils.nc'):
    """Load all Fourier modes of the vector potential from an Anvac NetCDF file."""

    from os import path, strerror
    from errno import ENOENT
    import netCDF4
    from numpy import array, empty, asarray

    if not path.exists(field_file):
        raise FileNotFoundError(ENOENT, strerror(ENOENT), field_file)

    rootgrp = netCDF4.Dataset(field_file, 'r')
    try:
        grid = grid_t()
        grid.R = array(rootgrp['R'])
        grid.Z = array(rootgrp['Z'])
        grid.nR = grid.R.size
        grid.nZ = grid.Z.size
        grid.R_min = float(grid.R[0])
        grid.R_max = float(grid.R[-1])
        grid.Z_min = float(grid.Z[0])
        grid.Z_max = float(grid.Z[-1])

        ncoil = len(rootgrp.dimensions['coil_number'])
        nmode = len(rootgrp.dimensions['ntor'])
        mode_numbers = asarray(range(nmode), dtype=int)

        def _allocate():
            return empty((nmode, ncoil, grid.nR, grid.nZ), dtype=complex)

        AnR = _allocate()
        Anphi = _allocate()
        AnZ = _allocate()
        dAnphi_dR = _allocate()
        dAnphi_dZ = _allocate()

        for n in range(nmode):
            for kcoil in range(ncoil):
                AnR[n, kcoil].real = rootgrp['AnR_real'][kcoil, :, :, n].T
                AnR[n, kcoil].imag = rootgrp['AnR_imag'][kcoil, :, :, n].T
                Anphi[n, kcoil].real = rootgrp['Anphi_real'][kcoil, :, :, n].T
                Anphi[n, kcoil].imag = rootgrp['Anphi_imag'][kcoil, :, :, n].T
                AnZ[n, kcoil].real = rootgrp['AnZ_real'][kcoil, :, :, n].T
                AnZ[n, kcoil].imag = rootgrp['AnZ_imag'][kcoil, :, :, n].T
                dAnphi_dR[n, kcoil].real = rootgrp['dAnphi_dR_real'][kcoil, :, :, n].T
                dAnphi_dR[n, kcoil].imag = rootgrp['dAnphi_dR_imag'][kcoil, :, :, n].T
                dAnphi_dZ[n, kcoil].real = rootgrp['dAnphi_dZ_real'][kcoil, :, :, n].T
                dAnphi_dZ[n, kcoil].imag = rootgrp['dAnphi_dZ_imag'][kcoil, :, :, n].T
    finally:
        rootgrp.close()

    return grid, mode_numbers, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ


def reconstruct_field_from_modes(
    BnR_modes,
    Bnphi_modes,
    BnZ_modes,
    mode_numbers,
    R_grid,
    Z_grid,
    x,
    y,
    z,
):
    """Reconstruct B(x, y, z) from Fourier modes.

    Parameters
    ----------
    BnR_modes, Bnphi_modes, BnZ_modes : ndarray
        Arrays of shape (nmodes, nR, nZ) or (nmodes, ...) containing coil-summed
        Fourier amplitudes.
    mode_numbers : array_like
        Sequence of non-negative toroidal mode numbers corresponding to modes.
    R_grid, Z_grid : ndarray
        1-D arrays describing the cylindrical grid.
    x, y, z : array_like
        Cartesian coordinates (cm) of evaluation points.

    Returns
    -------
    Bx, By, Bz : ndarray
        Real-valued Cartesian magnetic field components at the requested points.
    """

    from numpy import hypot, arctan2, cos, sin, exp, real, asarray, empty_like

    x = asarray(x, dtype=float)
    y = asarray(y, dtype=float)
    z = asarray(z, dtype=float)
    if x.shape != y.shape or x.shape != z.shape:
        raise ValueError("x, y, z must have identical shapes")

    flat_shape = x.shape
    x_flat = x.ravel()
    y_flat = y.ravel()
    z_flat = z.ravel()

    nmodes = len(mode_numbers)

    def _bilinear(field, R_val, Z_val):
        if R_val < R_grid[0] or R_val > R_grid[-1] or Z_val < Z_grid[0] or Z_val > Z_grid[-1]:
            return complex('nan')
        i1 = R_grid.searchsorted(R_val)
        if i1 == 0:
            i0 = i1 = 0
            t = 0.0
        elif i1 >= R_grid.size:
            i0 = i1 = R_grid.size - 1
            t = 0.0
        else:
            i0 = i1 - 1
            r0, r1 = R_grid[i0], R_grid[i1]
            t = 0.0 if r1 == r0 else (R_val - r0) / (r1 - r0)

        j1 = Z_grid.searchsorted(Z_val)
        if j1 == 0:
            j0 = j1 = 0
            u = 0.0
        elif j1 >= Z_grid.size:
            j0 = j1 = Z_grid.size - 1
            u = 0.0
        else:
            j0 = j1 - 1
            z0, z1 = Z_grid[j0], Z_grid[j1]
            u = 0.0 if z1 == z0 else (Z_val - z0) / (z1 - z0)

        f00 = field[i0, j0]
        f01 = field[i0, j1]
        f10 = field[i1, j0]
        f11 = field[i1, j1]
        return ((1.0 - t) * (1.0 - u) * f00 +
                (1.0 - t) * u * f01 +
                t * (1.0 - u) * f10 +
                t * u * f11)

    Bx = empty_like(x_flat)
    By = empty_like(y_flat)
    Bz = empty_like(z_flat)

    for idx, (x_val, y_val, z_val) in enumerate(zip(x_flat, y_flat, z_flat)):
        R_val = hypot(x_val, y_val)
        phi_val = arctan2(y_val, x_val)
        BR_total = 0.0 + 0.0j
        Bphi_total = 0.0 + 0.0j
        BZ_total = 0.0 + 0.0j

        for mode_idx, ntor in enumerate(mode_numbers):
            BR_mode = _bilinear(BnR_modes[mode_idx], R_val, z_val)
            Bphi_mode = _bilinear(Bnphi_modes[mode_idx], R_val, z_val)
            BZ_mode = _bilinear(BnZ_modes[mode_idx], R_val, z_val)

            if ntor == 0:
                BR_total += BR_mode
                Bphi_total += Bphi_mode
                BZ_total += BZ_mode
            else:
                factor = exp(1j * ntor * phi_val)
                BR_total += 2.0 * real(BR_mode * factor)
                Bphi_total += 2.0 * real(Bphi_mode * factor)
                BZ_total += 2.0 * real(BZ_mode * factor)

        cosphi = cos(phi_val)
        sinphi = sin(phi_val)
        Bx[idx] = (BR_total * cosphi - Bphi_total * sinphi).real
        By[idx] = (BR_total * sinphi + Bphi_total * cosphi).real
        Bz[idx] = BZ_total.real

    return Bx.reshape(flat_shape), By.reshape(flat_shape), Bz.reshape(flat_shape)


def gauged_Anvac_from_Bnvac(grid, BnR, BnZ, ntor=2):
    from numpy import newaxis

    gauged_AnR = 1j / ntor * grid.R[newaxis, :, newaxis] * BnZ
    gauged_AnZ = -1j / ntor * grid.R[newaxis, :, newaxis] * BnR
    return gauged_AnR, gauged_AnZ


def gauge_Anvac(grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=2):
    from numpy import newaxis

    if ntor == 0:
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
    }
    for kcoil in range(ncoil):
        spl['AnR_Re'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnR[kcoil].real, kx=5, ky=5)
        spl['AnR_Im'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnR[kcoil].imag, kx=5, ky=5)
        spl['AnZ_Re'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnZ[kcoil].real, kx=5, ky=5)
        spl['AnZ_Im'][kcoil] = RectBivariateSpline(grid.R, grid.Z, gauged_AnZ[kcoil].imag, kx=5, ky=5)

    # For ntor=0, also spline Aphi since we need it for B = curl(A)
    if ntor == 0 and Anphi is not None:
        spl['Anphi_Re'] = [None] * ncoil
        spl['Anphi_Im'] = [None] * ncoil
        for kcoil in range(ncoil):
            spl['Anphi_Re'][kcoil] = RectBivariateSpline(grid.R, grid.Z, Anphi[kcoil].real, kx=5, ky=5)
            spl['Anphi_Im'][kcoil] = RectBivariateSpline(grid.R, grid.Z, Anphi[kcoil].imag, kx=5, ky=5)

    return spl


def field_divfree(spl, R, Z, ntor=2):
    """Evaluate vector potential splines and return Fourier amplitude of magnetic field.

    Keyword arguments:
    spl -- dict of splines from spline_gauged_Anvac
    R -- array of radii in cylindrical coordinates
    Z -- array of altitudes in cylindrical coordinates

    Return values:
    BnR -- R component of the magnetic field's Fourier mode
    Bnphi -- phi component of the magnetic field's Fourier mode
    BnZ -- Z component of the magnetic field's Fourier mode

    The returned arrays have shape (ncoil, nR, nZ), where nR and nZ
    are the number of elements in the arrays R and Z, respectively.
    ncoil is the number of coils taken from the input splines.
    """

    from numpy import atleast_1d, empty, newaxis, squeeze

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
        # BR = -dAphi/dZ
        # Bphi = dAR/dZ - dAZ/dR
        # BZ = dAphi/dR + Aphi/R
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
        BnR = 1j * ntor * AnZ / R[newaxis, :, newaxis]
        Bnphi = dAnR_dZ - dAnZ_dR
        BnZ = -1j * ntor * AnR / R[newaxis, :, newaxis]
    return squeeze(BnR), squeeze(Bnphi), squeeze(BnZ)
