"""
libneo
-----------

Python interface and helpers for libneo.

Notes
- Keep top-level imports lightweight: avoid importing submodules with heavy
  optional dependencies at import time to prevent hard failures on systems
  without those deps (e.g. CI without MARS/OMFIT).
- Expose common submodules as attributes for backwards compatibility, e.g.
  ``from libneo import eqdsk``.
"""

# Re-export common submodules with graceful degradation on missing deps.

# Safe, pure-Python or compiled modules expected to exist after build
try:
    from . import eqdsk as eqdsk
except Exception:
    eqdsk = None

try:
    from . import eqdsk_base as eqdsk_base
except Exception:
    eqdsk_base = None

try:
    from . import flux_converter as flux_converter
except Exception:
    flux_converter = None

try:
    from . import interpolate as interpolate
except Exception:
    interpolate = None

try:
    from .semi_periodic_fourier_spline import (
        SemiPeriodicFourierSpline,
        SemiPeriodicFourierSpline_ExceedingNyquistLimit,
    )
except Exception:
    SemiPeriodicFourierSpline = None
    SemiPeriodicFourierSpline_ExceedingNyquistLimit = None

try:
    from .fourier_utils import (
        fourier_coefs_half,
        fourier_coefs_full,
        get_half_fft,
    )
except Exception:
    fourier_coefs_half = None
    fourier_coefs_full = None
    get_half_fft = None

# Optional/heavy dependencies: import lazily if available
try:
    from . import boozer as boozer
except Exception:
    boozer = None

try:
    from .coordinate_converter import (
        StorGeom2MarsCoords,
        order_monotonically,
    )
except Exception:
    StorGeom2MarsCoords = None
    order_monotonically = None

try:
    from .getHeaderDataVMEC import *  # small helper
except Exception:
    pass

try:
    from .convert_neo_in_to_fortran_namelist import *
except Exception:
    pass
