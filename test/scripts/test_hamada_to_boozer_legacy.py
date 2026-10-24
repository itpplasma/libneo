"""Legacy hamada_to_boozer.x on the public EFIT file test/resources/input_efit_file.dat.

The field-line integration of the converter must give the same iota as
efit_to_boozer.x on the same s grid, and an axisymmetric (n = 0) Hamada
spectrum with only m = 0 must convert to the same m = 0 Boozer amplitude with
all other modes zero. A shape mismatch in rhs_converter once returned
iota = 2 pi/h and NaN on every surface.
"""
import math
import pathlib
import shutil
import subprocess
import sys
import tempfile

efit_binary = pathlib.Path(sys.argv[1]).resolve()
hamada_binary = pathlib.Path(sys.argv[2]).resolve()
resources = pathlib.Path(__file__).resolve().parents[1] / "resources"
gfile = resources / "input_efit_file.dat"
nsurf, mpolmax, nrho, b0 = 20, 4, 21, 0.37


def blocks(path):
    """(s, iota, [(m, bmnc, bmns)]) per surface of a .bc file."""
    lines = path.read_text().splitlines()[6:]
    out = []
    while lines:
        profile = list(map(float, lines[2].split()))
        rows = []
        k = 4
        while k < len(lines) and not lines[k].lstrip().startswith("s "):
            values = lines[k].split()
            rows.append((int(values[0]), float(values[8]), float(values[9])))
            k += 1
        out.append((profile[0], profile[1], rows))
        lines = lines[k:]
    return out


with tempfile.TemporaryDirectory(dir=pathlib.Path.cwd()) as temporary:
    root = pathlib.Path(temporary)
    shutil.copy(gfile, root / gfile.name)
    shutil.copy(resources / "convexwall.dat", root / "convexwall.dat")
    (root / "field_divB0.inp").write_text(
        "0 ipert\n1 iequil\n0.00 ampl\n0 ntor\n0.99 cutoff\n4 icftype\n"
        f"'{gfile.name}' gfile\n'' pfile\n'convexwall.dat' convexfile\n'' fluxdatapath\n"
        "0 nwindow_r\n0 nwindow_z\n1 ieqfile\n")
    (root / "efit_to_boozer.inp").write_text(
        f"3600\n500\n1024\n10000\n{nsurf}\n24\n1e10\n2 construction (legacy)\n")
    subprocess.run([str(efit_binary)], cwd=root, check=True, capture_output=True)
    # MARS inputs: mpolmax, nrhopol, theta_H0, theta_B0, ntor; Hamada b_m = b0 delta_m0.
    (root / "mars_data_dimensions_and_extra.in").write_text(f"{mpolmax}\n{nrho}\n0.0\n0.0\n0\n")
    (root / "RMZM_F.OUT").write_text("1.0\n0.0 0.0 0.0 1.0\n")
    spectrum = [0.0] * (4 * mpolmax + 2)
    spectrum[mpolmax] = b0
    (root / "TORQUENTV.OUT").write_text("".join(f"{k / (nrho - 1):.16e}\n" for k in range(nrho)))
    (root / "DATATORQNTV.OUT").write_text((" ".join(map(repr, spectrum)) + "\n") * nrho)
    subprocess.run([str(hamada_binary)], cwd=root, check=True, capture_output=True)
    reference = blocks(root / "fromefit.bc")
    converted = blocks(root / "pert.bc")

assert len(converted) == len(reference) == nsurf
assert all(math.isfinite(c[1]) for c in converted), "non-finite iota"
worst_iota = max(abs(c[1] / r[1] - 1) for c, r in zip(converted, reference))
worst_mode = 0.0
for s, iota, rows in converted:
    for m, bc, bs in rows:
        expected = b0 if m == 0 else 0.0
        worst_mode = max(worst_mode, abs(bc - expected) / b0, abs(bs) / b0)
assert worst_iota < 1e-3, worst_iota
assert worst_mode < 1e-7, worst_mode
print(f"hamada_to_boozer legacy: max |iota/iota_efit - 1| = {worst_iota:.2e}, "
      f"m = 0 conversion error = {worst_mode:.2e}: PASS")
