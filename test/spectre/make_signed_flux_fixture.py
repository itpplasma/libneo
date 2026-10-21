"""Keep the existing fixture, replacing only the global flux metadata."""
import shutil
import sys

import h5py

shutil.copyfile(sys.argv[1], sys.argv[2])
with h5py.File(sys.argv[2], "r+") as handle:
    handle["input/physics/phiedge"][...] = -3.25
