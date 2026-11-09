"""Pinned test data: fetched once per machine, retried, checked by SHA-256.

fetch(url) returns a local copy of a file listed in PINNED. Copies live in
LIBNEO_TEST_DATA_CACHE (default: a per-user directory in the temp dir) and
are reused while their hash matches.
"""
from __future__ import annotations

import hashlib
import os
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

SIMSOPT = ("https://raw.githubusercontent.com/hiddenSymmetries/simsopt/"
           "912bdcc729d95433518c71df82b8cc45287908b3/tests/test_files/")
BOOZ_XFORM = ("https://raw.githubusercontent.com/hiddenSymmetries/booz_xform/"
              "21aff118735de4255641aa1b64fc7ab35e8c9767/tests/test_files/")
STELLOPT = "https://princetonuniversity.github.io/STELLOPT/examples/"

WOUT_QA_URL = SIMSOPT + "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
WOUT_QH_URL = SIMSOPT + "wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
WOUT_LSP_URL = SIMSOPT + "wout_LandremanSenguptaPlunk_section5p3_reference.nc"
WOUT_LI383_URL = BOOZ_XFORM + "wout_li383_1.4m.nc"
WOUT_CIRCULAR_TOKAMAK_URL = BOOZ_XFORM + "wout_circular_tokamak.nc"
WOUT_UPDOWN_TOKAMAK_URL = BOOZ_XFORM + "wout_up_down_asymmetric_tokamak.nc"
WOUT_NCSX_URL = STELLOPT + "wout_ncsx_c09r00_fixed.nc"
COILS_NCSX_URL = STELLOPT + "coils.c09r00"

PINNED = {
    WOUT_QA_URL: "333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771",
    WOUT_QH_URL: "249c397d67357949b2ab4690145224736d8388f5332f23c67c9cb0bf2aa9931e",
    WOUT_LSP_URL: "a5bc177285779cb11a9dc4c7d95d3392ffc295b44159e47b35f581587596bc98",
    WOUT_LI383_URL: "6508955581628d8a67992d340fe3a82c160c1a183925d8dfb1f955a4dd596a3f",
    WOUT_CIRCULAR_TOKAMAK_URL:
        "625c5de8e0eb6041b505ffee17e13d3c1b1837488a6c3d1c9b351ca06bc6749a",
    WOUT_UPDOWN_TOKAMAK_URL:
        "5fb6f1c27f84a1ba226a9509e17d383de60874175293c50fdcd89c51649fd1a9",
    WOUT_NCSX_URL: "1005d81729366118ed51926c51b99b08bea21cfa77bb6f2874774632fb0a1e99",
    COILS_NCSX_URL: "3c429da06f4c062887a497a16e2d2bd10f0ecb0b8858c252698631f3853da428",
}


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cache_dir() -> Path:
    default = Path(tempfile.gettempdir()) / f"libneo-test-data-{os.getuid()}"
    return Path(os.environ.get("LIBNEO_TEST_DATA_CACHE", default))


def download(url: str, path: Path, sha256: str, attempts: int = 5) -> Path:
    """Make path hold the file at url with the given SHA-256."""
    path = Path(path)
    if path.exists() and sha256sum(path) == sha256:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(f"{path.name}.{os.getpid()}.part")
    error = ""
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                part.write_bytes(response.read())
            actual = sha256sum(part)
            if actual == sha256:
                part.replace(path)
                return path
            error = f"SHA-256 {actual}, expected {sha256}"
        except (urllib.error.URLError, OSError) as exc:
            error = str(exc)
        if attempt < attempts:
            time.sleep(10 * attempt)
    part.unlink(missing_ok=True)
    raise RuntimeError(f"cannot fetch {url}: {error}")


def fetch(url: str) -> Path:
    """Cached, verified local copy of a pinned file."""
    sha = PINNED[url]
    return download(url, cache_dir() / sha[:16] / url.rsplit("/", 1)[-1], sha)
