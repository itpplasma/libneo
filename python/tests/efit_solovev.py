"""Analytic synthetic Solov'ev equilibrium shared by EFIT consumer tests.

The polynomial and independent contour quadrature are the same mathematical
oracle used by test/source/test_efit_sfl.f90. No machine data are included.
"""
import numpy as np


R0, B0, A, B, C, D, PSI_EDGE = 1.7, 2.0, 0.092, 0.079, 0.2, 0.01, 0.2
MU0 = 4e-7 * np.pi


def _psi(r, z):
    return A*r*r*z*z + B*(r*r-R0*R0)**2 + C*z*z + D*r*r*z


def _grad(r, z):
    return (2*A*r*z*z + 4*B*r*(r*r-R0*R0) + 2*D*r*z,
            2*A*r*r*z + 2*C*z + D*r*r)


def _axis():
    r, z = R0, 0.0
    for _ in range(12):
        hessian = [[2*A*z*z + 12*B*r*r - 4*B*R0*R0 + 2*D*z,
                    4*A*r*z + 2*D*r],
                   [4*A*r*z + 2*D*r, 2*A*r*r + 2*C]]
        dr, dz = np.linalg.solve(hessian, _grad(r, z))
        r, z = r-dr, z-dz
    return np.array([r, z])


def _contour_oracle(psi, psi_sign=1, f_sign=1):
    """Return q_CCW, I_CCW and G in SI by polar ray bisection."""
    angle = np.arange(1024)*2*np.pi/1024
    co, si = np.cos(angle), np.sin(angle)
    raxis, zaxis = _axis()
    low, high = np.zeros(angle.size), np.full(angle.size, 1.5)
    for _ in range(65):
        middle = (low+high)/2
        outside = _psi(raxis+middle*co, zaxis+middle*si) > psi
        high = np.where(outside, middle, high)
        low = np.where(outside, low, middle)
    rho = (low+high)/2
    r, z = raxis+rho*co, zaxis+rho*si
    pr, pz = _grad(r, z)
    radial_grad = pr*co+pz*si
    f = f_sign*np.sqrt((R0*B0)**2 - 4*C*(psi-PSI_EDGE))
    q = psi_sign*f*np.mean(rho/(r*radial_grad))
    current_i = psi_sign*np.mean((pr*pr+pz*pz)*rho/(r*radial_grad))
    return q, current_i, f


def _write_solovev(path, psi_sign=1, f_sign=1, pprime_scale=1):
    n = 65
    r, z = np.linspace(.9, 2.4, n), np.linspace(-.9, .9, n)
    raxis, zaxis = _axis()
    psi_axis = _psi(raxis, zaxis)
    flux = np.linspace(psi_axis, PSI_EDGE, n)
    q = [_contour_oracle(psi_axis + max(x, 1e-6)*(PSI_EDGE-psi_axis))[0]
         for x in np.linspace(0, 1, n)]
    rr, zz = np.meshgrid(r, z)
    pprime = -(8*B+2*A)/MU0
    header = [1.5, 1.8, R0, .9, 0, raxis, zaxis,
              psi_sign*psi_axis, psi_sign*PSI_EDGE, f_sign*B0,
              _contour_oracle(PSI_EDGE, psi_sign, f_sign)[1]*5e6, psi_sign*psi_axis, 0, raxis, 0,
              zaxis, 0, psi_sign*PSI_EDGE, 0, 0]
    arrays = [header,
              f_sign*np.sqrt((R0*B0)**2-4*C*(flux-PSI_EDGE)),
              pprime_scale*pprime*(flux-PSI_EDGE),
              np.full(n, -2*C*psi_sign),
              np.full(n, pprime_scale*pprime*psi_sign),
              psi_sign*_psi(rr, zz).ravel(), q]
    with path.open("w") as stream:
        stream.write(f"{'analytic Solovev, synthetic':48s}{0:4d}{n:4d}{n:4d}\n")
        for values in arrays:
            for i in range(0, len(values), 5):
                stream.write("".join(f"{v:16.9E}" for v in values[i:i+5])+"\n")
        stream.write(f"{1:5d}{0:5d}\n{raxis:16.9E}{zaxis:16.9E}\n")
