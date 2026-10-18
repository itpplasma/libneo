import os
import tempfile
from urllib.request import urlopen

import pytest

import numpy as np
from libneo.coils import CoilsFile, Filament


STELLOPT_COILS_URL = (
    "https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00"
)


def _download_text(url: str) -> str:
    with urlopen(url, timeout=30) as resp:
        return resp.read().decode("utf-8")


@pytest.mark.network
def test_parse_stellopt_example_roundtrip():
    text = _download_text(STELLOPT_COILS_URL)
    assert len(text) > 0

    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "coils.c09r00")
        with open(src, "w", encoding="utf-8") as f:
            f.write(text)

        cf = CoilsFile.from_file(src)
        assert cf.periods in (None, 3)
        assert len(cf.filaments) >= 1
        f0 = cf.filaments[0]
        assert f0.coords.shape[1] == 3
        assert f0.coords.shape[0] == f0.current.shape[0]
        # Some current values should be non-zero
        assert (f0.current != 0.0).any()

        # Map to Mgrid metadata (averages currents per filament)
        groups, currents = cf.to_mgrid_metadata()
        assert len(groups) == len(currents) == len(cf.filaments)
        assert all(isinstance(g, str) for g in groups)
        assert all(isinstance(c, float) for c in currents)

        # Round-trip write/read
        dst = os.path.join(td, "roundtrip.coils")
        cf.write(dst)
        cf2 = CoilsFile.from_file(dst)

        assert cf2.periods == cf.periods
        assert len(cf2.filaments) == len(cf.filaments)
        import numpy as np
        for a, b in zip(cf.filaments, cf2.filaments):
            assert np.allclose(a.coords, b.coords)
            assert np.allclose(a.current, b.current)

        # Test write_simple format
        simple_file = os.path.join(td, "coils.simple")
        cf.write_simple(simple_file)

        # Verify simple format structure
        with open(simple_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            n_points = int(lines[0].strip())
            assert n_points == sum(len(fil.coords) + 1 for fil in cf.filaments)  # + closing point
            assert len(lines) == n_points + 1  # header + data lines

            # Check first data line has 4 numbers
            parts = lines[1].split()
            assert len(parts) == 4
            for part in parts:
                float(part)  # should not raise


def _square_loop(z, side=1.0):
    h = side / 2.0
    pts = [(-h, -h, z), (h, -h, z), (h, h, z), (-h, h, z)]
    return Filament(coords=np.array(pts), current=np.full(4, 1000.0))


def _read_simple(path):
    with open(path, encoding="utf-8") as f:
        n = int(f.readline())
        data = np.array([[float(v) for v in f.readline().split()] for _ in range(n)])
    return data


def _chain_field(data, x):
    """Biot-Savart of a polygonal chain; segment i carries current of point i."""
    b = np.zeros(3)
    for p, q, cur in zip(data[:-1, :3], data[1:, :3], data[:-1, 3]):
        r1, r2 = x - p, x - q
        c = np.cross(q - p, r1)
        if cur == 0.0 or np.dot(c, c) == 0.0:
            continue
        b += cur * c / np.dot(c, c) * (
            np.dot(q - p, r1) / np.linalg.norm(r1)
            - np.dot(q - p, r2) / np.linalg.norm(r2)
        )
    return 1e-7 * b  # SI, mu0/4pi = 1e-7


def _square_loop_axis_field(current, a, z):
    """Analytic on-axis Bz [T] of a square loop of side a at distance z."""
    mu0 = 4e-7 * np.pi
    return mu0 * current * a**2 / (2.0 * np.pi * (z**2 + a**2 / 4.0)
                                   * np.sqrt(z**2 + a**2 / 2.0))


def test_write_simple_closes_coils_and_separates_with_zero_current(tmp_path):
    coils = CoilsFile(periods=None, mirror=None,
                      filaments=[_square_loop(0.0), _square_loop(5.0)])
    out = tmp_path / "coils.simple"
    coils.write_simple(str(out))
    data = _read_simple(out)

    assert data.shape == (10, 4)
    assert np.array_equal(data[4, :3], data[0, :3]) and data[4, 3] == 0.0
    assert np.array_equal(data[9, :3], data[5, :3]) and data[9, 3] == 0.0
    assert np.count_nonzero(data[:, 3] == 0.0) == 2
    seg_current = data[:-1, 3]
    assert np.all(seg_current[[4]] == 0.0)  # segment between the loops
    assert np.all(seg_current[[0, 1, 2, 3, 5, 6, 7, 8]] == 1000.0)


def test_write_simple_field_matches_analytic_square_loop(tmp_path):
    a, d, cur = 1.0, 5.0, 1000.0
    coils = CoilsFile(periods=None, mirror=None,
                      filaments=[_square_loop(0.0, a), _square_loop(d, a)])
    out = tmp_path / "coils.simple"
    coils.write_simple(str(out))
    data = _read_simple(out)

    x = np.array([0.0, 0.0, 0.0])
    expected = _square_loop_axis_field(cur, a, 0.0) + _square_loop_axis_field(cur, a, d)
    bz = _chain_field(data, x)[2]
    assert bz == pytest.approx(expected, rel=1e-10)
