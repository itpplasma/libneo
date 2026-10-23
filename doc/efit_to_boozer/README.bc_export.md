# Checked G-EQDSK to Boozer files

`libneo.convert_eqdsk_to_bc` writes an asymmetric, left-handed `.bc` file
from the native contour construction. It rejects failed equilibrium checks
and nonfinite data before replacing the destination file.

```python
from libneo.eqdsk_to_boozer import convert_eqdsk_to_bc

convert_eqdsk_to_bc("equilibrium.g", "equilibrium.bc")
```

The optional keywords are `nsurf=64`, `mpol=32`, `n_theta=1024`, `n_q=128`,
`route=1`, `psimax=None`, and `exp_bpol=-1`. The surfaces are midpoints of
uniform intervals in normalized toroidal flux. Routes 1 and 2 construct
Boozer angles. `psimax` selects an optional axis-relative poloidal flux
boundary in CGS. `exp_bpol` selects automatic flux-unit inference (-1),
Wb/rad (0), or Wb (1). Quality overrides are unavailable in this exporter.

Geometry and magnetic-field spectra, currents, pressure derivative, volume
derivative, axis and minor radius all come from the same checked native
transform. The file uses SI units and the existing LHS current and signed
volume conventions. The native angular displacement remains dimensionless.

The exporter requests `precision=16` from the existing Boozer writers.
Their default output format remains unchanged. They also accept this
keyword directly, including `BoozerFile.write(path, precision=16)`.
