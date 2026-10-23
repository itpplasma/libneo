"""Legacy field-line turn detection on a low-q, strongly shaped equilibrium.

test/resources/freegs_mast_129.geqdsk is public and synthetic: FreeGS 5b41fe5
example 03 (MAST), 129 x 129, rtol 1e-8, profile columns written on the
G-EQDSK grid psi_N = k/(nx - 1) (sha256 58d225c3...). With q near 0.5 a poloidal
turn takes about five scan steps, and plain Newton from the end of the step
that crossed Z = zaxis could converge to the inboard crossing, giving q of the
wrong size or sign on some start radii. The file's own qpsi is an independent
oracle (it agrees with contour q to 8e-4 inside 0.99 psi_sep, and linear
interpolation adds up to 2e-3 for psi_N <= 0.95). The defect gave errors above 1.
"""
import pathlib
import shutil
import subprocess
import sys
import tempfile

binary = pathlib.Path(sys.argv[1]).resolve()
resources = pathlib.Path(__file__).resolve().parents[1] / "resources"
gfile = resources / "freegs_mast_129.geqdsk"

tokens = gfile.read_text().split("\n", 1)[1].replace("E-", "E@").replace("-", " -").replace("E@", "E-").split()
values = list(map(float, tokens))
nw = int(gfile.read_text().split("\n", 1)[0].split()[-2])
simag, sibry = values[7], values[8]
qpsi = values[20 + 4 * nw + nw * nw:20 + 5 * nw + nw * nw]

with tempfile.TemporaryDirectory(dir=pathlib.Path.cwd()) as temporary:
    root = pathlib.Path(temporary)
    shutil.copy(gfile, root / gfile.name)
    shutil.copy(resources / "convexwall.dat", root / "convexwall.dat")
    (root / "field_divB0.inp").write_text(
        "0 ipert\n1 iequil\n0.00 ampl\n0 ntor\n0.99 cutoff\n4 icftype\n"
        f"'{gfile.name}' gfile\n'' pfile\n'convexwall.dat' convexfile\n'' fluxdatapath\n"
        "0 nwindow_r\n0 nwindow_z\n1 ieqfile\n")
    psimax = 0.99 * (sibry - simag) * 1e8
    (root / "efit_to_boozer.inp").write_text(
        f"3600\n500\n500\n10000\n100\n40\n{psimax:.17g}\n")
    subprocess.run([str(binary)], cwd=root, check=True, capture_output=True)
    rows = [list(map(float, line.split()))
            for line in (root / "flux_functions.dat").read_text().splitlines()[1:]]

worst = 0.0
for row in rows:
    if row[3] * 0.99 > 0.95:  # linear interpolation of qpsi is too coarse at the edge
        continue
    psin = row[3] * 0.99 * (nw - 1)  # psi_pol is normalized to psimax
    k = min(int(psin), nw - 2)
    q_file = abs(qpsi[k] + (psin - k) * (qpsi[k + 1] - qpsi[k]))
    worst = max(worst, abs(row[2] / q_file - 1.0))
assert worst < 1e-2, worst
print(f"Low-q field-line turn detection: max |q/qpsi - 1| = {worst:.2e}: PASS")
