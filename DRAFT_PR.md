# Draft PR: Add geoflux coordinate support for GEQDSK files

## Summary
- add `geoflux_coordinates` module that reads a GEQDSK equilibrium, builds the toroidal flux mapping, and exposes `geoflux↔cyl/cart` transforms reusing the simple VMEC-style interface
- wire the new module into the core library build and provide a lightweight regression exercising round-trip conversions against the existing MAST GEQDSK sample
- extend the CMake test harness with `test_geoflux.x`

## Testing
- `cmake --build build -j` *(fails: CMake configuration requires Python development headers with NumPy; environment lacks `Python3_NumPy`)*
