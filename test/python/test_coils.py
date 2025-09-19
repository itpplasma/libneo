import os
import tempfile
from urllib.request import urlopen

import pytest

from libneo.coils import CoilsFile


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
