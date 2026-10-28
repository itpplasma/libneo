# libneo Python package

The Python package lives in `python/libneo` (with `python/efit_to_boozer`). It is
built with scikit-build-core; the Fortran wrappers are built only if `numpy`
and `f90wrap` are available. From the repository root, inside a virtual
environment:

    pip install -e .                 # base package
    pip install -e ".[chartmap]"     # chartmap generators (map2disc, shapely, trimesh, ...)

Command-line entry points defined in `pyproject.toml`: `libneo-write-chartmap`,
`libneo-booz-xform`, `nemec2vmec` and `neo2bc2chartmap`.

## Chartmap generators

These write `*.chartmap.nc` volumes that can be read by the Fortran module
`libneo_coordinates`.

For map2disc-generated chartmaps, `zeta` is the cylindrical toroidal angle
(`zeta_convention=cyl`) and `rho` is a geometric disk radius
(`rho_convention=unknown`).

Python API:

    from libneo.chartmap import write_chartmap_from_vmec_boundary
    write_chartmap_from_vmec_boundary("wout.nc", "wout.chartmap.nc", nrho=33, ntheta=65, nzeta=33)

Command line:

    libneo-write-chartmap from-vmec wout.nc wout.chartmap.nc --nrho 33 --ntheta 65 --nzeta 33

### VMEC with extension beyond the last closed flux surface

    libneo-write-chartmap from-vmec-extended wout.nc wout.chartmap.nc \
        --nrho-vmec 25 --boundary-offset 0.15

Inside the LCFS the exact VMEC coordinates are used. Outside, coordinate lines
are continued by cubic Hermite interpolation along their tangents, which gives
C1 continuity at the boundary, as needed for symplectic integrators. No map2disc
dependency is required.

### VMEC extended to an STL wall

    libneo-write-chartmap from-vmec-to-wall wout.nc wall.stl wout.chartmap.nc \
        --rho-lcfs 0.8 --nrho 33 --ntheta 65 --nzeta 33

VMEC coordinates are used inside the LCFS and linearly interpolated to a wall
contour extracted from the STL file. The radial coordinate `rho` runs from 0 to
1; `rho <= rho_lcfs` corresponds to VMEC surfaces with `s = (rho/rho_lcfs)^2`,
and `rho = 1` lies on the wall. Poloidal points on the LCFS and the wall are
matched by starting both curves at the outboard midplane and interpolating by
fractional arc position.

STL files carry no unit information; meters are assumed. For CAD files in
millimeters pass `--stl-units mm` (or `--stl-scale 1e-3`).

## VMEC angle conventions

VMEC documentation usually writes Fourier factors as `exp(i*(m*theta - n*zeta))`.
In VMEC `wout_*.nc` files

- `xm` stores the poloidal mode numbers `m`;
- `xn` stores signed toroidal mode numbers including the field periods,
  `xn = n*nfp`;
- `zeta` (or `varphi`) is the geometric toroidal angle in radians. Some tools
  reduce it modulo one field period, which does not change the Fourier phase
  since `xn` already contains `nfp`.

Real-space evaluation with `wout` coefficients therefore uses
`cos(m*theta - xn*zeta)` and `sin(m*theta - xn*zeta)`. The magnetic axis uses
the same phase, so `Z_axis(zeta) = sum_n zaxis_cs(n) * sin(-n*nfp*zeta)`.

This convention agrees with [STELLOPT `read_wout_mod.f90`][stellopt-wout] and
[DESC `vmec.py`][desc-vmec].

[stellopt-wout]: https://github.com/PrincetonUniversity/STELLOPT/blob/develop/LIBSTELL/Sources/Modules/read_wout_mod.f90
[desc-vmec]: https://github.com/PlasmaControl/DESC/blob/master/desc/vmec.py

## NEMEC to VMEC conversion

`nemec2vmec` converts text NEMEC `wout.*` files to VMEC-style NetCDF `wout.nc`:

    nemec2vmec wout.nemec wout.nc

The NEMEC output manual distributed with the Strumberger/Tichmann data
([nemec_output_user_manual.pdf](efit_to_boozer/nemec_output_user_manual.pdf))
uses the same geometric phase, `cos(m*theta - n*zeta*nfp)` and
`sin(m*theta - n*zeta*nfp)` for `R` and `Z`. Differences handled by the
converter:

- NEMEC stores `hphip = (1/(2*pi)) dPhi/ds`, while VMEC `phipf` and `phips`
  store `dPhi/ds`. For AUG-like equilibria where `Phi` decreases from axis to
  edge, `nemec2vmec` writes `-2*pi*hphip`.
- `signgs` is written as `-1`, matching the negative VMEC Jacobian convention.

## Coil format converter

Convert STELLOPT coil files to the simple Biot-Savart format read by the
Fortran module `neo_biotsavart` and by SIMPLE:

    python -m libneo.convert_coils_to_simple input.coils output.simple

## ASCOT5 helpers

- `libneo.ascot5.field_from_vmec` samples a VMEC equilibrium through the
  f2py-generated `_magfie` wrappers and returns an ASCOT5 `B_3DS` field (CGS
  internally, converted to SI on export).
- `libneo.ascot5.field_from_mgrid` converts an mgrid NetCDF file into the same
  `B3DSField` dataclass.
- `libneo.ascot5.write_b3ds_hdf5` writes the field to HDF5 with the layout
  expected by ASCOT5.

See [ascot5.md](ascot5.md) for the workflow and troubleshooting.
