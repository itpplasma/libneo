#!/usr/bin/env python3
"""Generate a tilted circular coil geometry for validation workflows."""

import numpy as np
from numpy import pi, cos, sin
import sys

def _plane_basis(center_phi, tilt_theta_deg, tilt_psi_deg):
    """Return orthonormal basis (u, v, n) for the coil plane."""
    theta = np.deg2rad(tilt_theta_deg)
    psi = np.deg2rad(tilt_psi_deg)
    phi0 = center_phi

    e_R = np.array([cos(phi0), sin(phi0), 0.0])
    e_phi = np.array([-sin(phi0), cos(phi0), 0.0])
    e_Z = np.array([0.0, 0.0, 1.0])

    n_vec = sin(theta) * cos(psi) * e_R + sin(theta) * sin(psi) * e_phi + cos(theta) * e_Z
    n_vec = n_vec / np.linalg.norm(n_vec)

    # Choose first in-plane direction to avoid degeneracy
    candidate = np.cross(n_vec, e_Z)
    if np.linalg.norm(candidate) < 1e-12:
        candidate = np.cross(n_vec, e_R)
    u_vec = candidate / np.linalg.norm(candidate)
    v_vec = np.cross(n_vec, u_vec)

    return u_vec, v_vec, n_vec


def create_tilted_coil(radius=2.0, center=(2.0, 0.35, 0.8),
                       tilt_theta=30.0, tilt_psi=35.0, num_segments=129):
    """Create a circular coil with general tilt relative to cylindrical axes.

    Parameters
    ----------
    radius : float
        Coil radius in meters.
    center : tuple
        Cylindrical coordinates (R, phi, Z) of the coil center.
    tilt_theta : float
        Polar tilt (deg) of coil normal relative to +Z.
    tilt_psi : float
        Azimuthal tilt (deg) of coil normal measured from +R toward +phi.
    num_segments : int
        Number of coil segments to discretize circle.
    """
    R0, phi0, Z0 = center
    u_vec, v_vec, _ = _plane_basis(phi0, tilt_theta, tilt_psi)

    # Coil center in Cartesian
    center_xyz = np.array([R0 * cos(phi0), R0 * sin(phi0), Z0])

    # Parametric circle
    t = np.linspace(0.0, 2.0 * pi, num_segments, endpoint=False)
    offsets = radius * (np.outer(u_vec, np.cos(t)) + np.outer(v_vec, np.sin(t)))
    points = center_xyz[:, None] + offsets
    return points

if __name__ == '__main__':
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'tilted_coil.dat'

    xyz = create_tilted_coil(
        radius=2.0,
        center=(2.0, 0.35, 0.8),
        tilt_theta=30.0,
        tilt_psi=35.0,
        num_segments=129,
    )

    nseg = xyz.shape[1]
    with open(output_file, 'w') as f:
        f.write(f"{nseg}\n")
        for i in range(nseg):
            # Current is 1.0 abampere (= 10 A) per segment
            f.write(f"{xyz[0,i]:23.16e} {xyz[1,i]:23.16e} {xyz[2,i]:23.16e} 1.0\n")

    print(f"Wrote {nseg}-segment tilted coil to {output_file}")
