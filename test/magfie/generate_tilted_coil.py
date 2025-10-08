#!/usr/bin/env python3
"""Generate a tilted circular coil for ntor=0 testing."""

import numpy as np
from numpy import pi, cos, sin
import sys

def create_tilted_coil(radius=2.0, center=(2.0, 0.0, 0.8),
                       tilt_angle=30.0, num_segments=129):
    """Create a circular coil tilted in the R-Z plane."""
    theta_rad = np.deg2rad(tilt_angle)
    R0, phi0, Z0 = center

    # Parametric angle around coil
    t = np.linspace(0, 2*pi, num_segments, endpoint=False)

    # Local coil coordinates (circular in tilted plane)
    x_local = radius * cos(t)
    y_local = radius * sin(t)

    # Coil normal is tilted in R-Z plane
    # Coil plane passes through (R0, Z0)
    x = R0 + x_local * cos(theta_rad)
    y = y_local
    z = Z0 + x_local * sin(theta_rad)

    # Convert to Cartesian
    X = x * cos(phi0) - y * sin(phi0)
    Y = x * sin(phi0) + y * cos(phi0)
    Z = z

    return np.vstack([X, Y, Z])

if __name__ == '__main__':
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'tilted_coil.dat'

    xyz = create_tilted_coil(radius=2.0, center=(2.0, 0.0, 0.8),
                             tilt_angle=30.0, num_segments=129)

    nseg = xyz.shape[1]
    with open(output_file, 'w') as f:
        f.write(f"{nseg}\n")
        for i in range(nseg):
            # Current is 1.0 abampere (= 10 A) per segment
            f.write(f"{xyz[0,i]:23.16e} {xyz[1,i]:23.16e} {xyz[2,i]:23.16e} 1.0\n")

    print(f"Wrote {nseg}-segment tilted coil to {output_file}")
