"""
libneo
-----------

This is the package containing the modules offered by the libneo
Python interface. Currently this is mainly the magfie module with support for
VMEC equilibrium fields. The magfie module is directly generated via f2py
on package install.
"""

from .boozer import *
from .eqdsk import *
from .eqdsk_base import *
from .convert_neo_in_to_fortran_namelist import *
from .getHeaderDataVMEC import *
from .flux_converter import *
from .biotsavart_fourier import *
from .semi_periodic_fourier_spline import SemiPeriodicFourierSpline
from .semi_periodic_fourier_spline import SemiPeriodicFourierSpline_ExceedingNyquistLimit
from .fourier_utils import fourier_coefs_half, fourier_coefs_full
from .fourier_utils import get_half_fft

from .coordinate_converter import StorGeom2MarsCoords
from .coordinate_converter import order_monotonically

# VMEC utilities
from .vmec import VMECGeometry, vmec_to_cylindrical

# Chartmap coordinate utilities
from .chartmap import chartmap_to_rz, vmec_to_chartmap_coords

# ASCOT5 utilities
from .ascot5 import field_from_vmec, field_from_mgrid, write_b3ds_hdf5
