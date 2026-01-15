EQDSK Equilibrium Files
=======================

The ``eqdsk_file`` class provides tools for reading G-EQDSK equilibrium files
and converting between cylindrical (R, Z) coordinates and flux coordinates.

Reading EQDSK Files
-------------------

.. code-block:: python

   from libneo import eqdsk_file

   eq = eqdsk_file('equilibrium.geqdsk')

   # Access grid data
   print(f"R range: {eq.R.min():.2f} - {eq.R.max():.2f} m")
   print(f"Z range: {eq.Z.min():.2f} - {eq.Z.max():.2f} m")
   print(f"Magnetic axis: R={eq.Rpsi0:.3f} m, Z={eq.Zpsi0:.3f} m")

Coordinate Conversion
---------------------

Converting (R, Z) to Flux Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most common use case is converting arbitrary (R, Z) measurement points
to normalized poloidal flux coordinates for radial profile analysis:

.. code-block:: python

   import numpy as np
   from libneo import eqdsk_file

   eq = eqdsk_file('equilibrium.geqdsk')

   # Single point
   R, Z = 1.8, 0.1
   s_pol, theta = eq.rz_to_flux_coords(R, Z)
   print(f"s_pol = {s_pol:.3f}, theta = {theta:.3f} rad")

   # Array of measurement points
   R_data = np.array([1.7, 1.8, 1.9, 2.0])
   Z_data = np.array([0.0, 0.1, -0.1, 0.0])
   s_pol, theta = eq.rz_to_flux_coords(R_data, Z_data)

Interpolating Poloidal Flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get the raw poloidal flux values at arbitrary (R, Z) points:

.. code-block:: python

   # Scalar
   psi = eq.psi_at_rz(1.8, 0.1)

   # Vectorized
   psi = eq.psi_at_rz(R_data, Z_data)

   # On a 2D grid (for contour plots)
   R_grid = np.linspace(1.5, 2.2, 50)
   Z_grid = np.linspace(-0.5, 0.5, 50)
   psi_2d = eq.psi_at_rz(R_grid, Z_grid, grid=True)

Plotting Scalar Data vs Flux Coordinate
---------------------------------------

A typical workflow for plotting diagnostic data against normalized flux:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from libneo import eqdsk_file

   # Load equilibrium
   eq = eqdsk_file('equilibrium.geqdsk')

   # Your measurement data at (R, Z) locations
   R_meas = np.array([1.65, 1.70, 1.75, 1.80, 1.85, 1.90])
   Z_meas = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
   temperature = np.array([3.2, 2.8, 2.3, 1.7, 1.1, 0.5])  # keV

   # Convert to flux coordinates
   s_pol, theta = eq.rz_to_flux_coords(R_meas, Z_meas)
   rho_pol = np.sqrt(s_pol)  # sqrt(s_pol) is common radial coordinate

   # Plot
   plt.figure()
   plt.plot(rho_pol, temperature, 'o-')
   plt.xlabel(r'$\\rho_{pol}$')
   plt.ylabel('Temperature [keV]')
   plt.title('Temperature Profile')
   plt.show()

Flux Coordinate Definitions
---------------------------

Normalized Poloidal Flux (s_pol)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The normalized poloidal flux is defined as:

.. math::

   s_{pol} = \\frac{\\psi - \\psi_{axis}}{\\psi_{edge} - \\psi_{axis}}

where:

- :math:`s_{pol} = 0` at the magnetic axis
- :math:`s_{pol} = 1` at the last closed flux surface (LCFS/separatrix)

Geometric Poloidal Angle (theta)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The geometric poloidal angle is defined as:

.. math::

   \\theta = \\arctan2(Z - Z_{axis}, R - R_{axis})

where:

- :math:`\\theta = 0` at the outboard midplane (low-field side)
- :math:`\\theta = \\pi/2` at the top
- :math:`\\theta = -\\pi/2` at the bottom
- :math:`\\theta = \\pm\\pi` at the inboard midplane (high-field side)

COCOS Conventions
-----------------

EQDSK files from different sources may use different coordinate conventions
(COCOS). The ``eqdsk_file`` class automatically detects sign inconsistencies
between the header values and the 2D psi array, ensuring that ``spol_at_rz``
always returns 0 at the axis and 1 at the LCFS regardless of the input file
convention.

API Reference
-------------

.. autoclass:: libneo.eqdsk.eqdsk_file
   :members: psi_at_rz, spol_at_rz, theta_geometric_at_rz, rz_to_flux_coords
   :show-inheritance:
