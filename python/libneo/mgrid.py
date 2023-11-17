"""
Handle mgrid data
"""
import textwrap
from netCDF4 import Dataset


class MgridFile:
    def __init__(self, filename=None):
        if filename:
            self.read(filename)

    def __str__(self):
        return textwrap.dedent(
            f"""
            mgrid
            ir = {self.ir}
            jz = {self.jz}
            kp = {self.kp}
            currents: {self.coil_current}
        """
        )

    def read(self, filename):
        with Dataset(filename, "r") as f:
            self.coil_group = f.variables["coil_group"][:]
            self.coil_current = f.variables["raw_coil_cur"][:]
            self.ir = f.dimensions["rad"].size
            self.jz = f.dimensions["zee"].size
            self.kp = f.dimensions["phi"].size

            self.nfp = f.variables["nfp"][:]
            self.rmin = f.variables["rmin"][:]
            self.zmin = f.variables["zmin"][:]
            self.rmax = f.variables["rmax"][:]
            self.zmax = f.variables["zmax"][:]
            self.mgrid_mode = f.variables["mgrid_mode"][:]

            self.br = []
            self.bp = []
            self.bz = []
            for kgroup, group in enumerate(self.coil_group):
                self.br.append(f.variables[f"br_{kgroup+1:03}"][:])
                self.bp.append(f.variables[f"bp_{kgroup+1:03}"][:])
                self.bz.append(f.variables[f"bz_{kgroup+1:03}"][:])

    def write(self, filename):
        with Dataset(filename, "w") as f:
            # Add dimensions
            f.createDimension("stringsize", 30)
            f.createDimension("external_coil_groups", self.coil_group.shape[0])
            f.createDimension("dim_00001", 1)
            f.createDimension("external_coils", self.coil_current.size)
            f.createDimension("rad", self.ir)
            f.createDimension("zee", self.jz)
            f.createDimension("phi", self.kp)

            # Add variables
            f.createVariable("ir", "i4", ())
            f.createVariable("jz", "i4", ())
            f.createVariable("kp", "i4", ())
            f.createVariable("nfp", "i4", ())
            f.createVariable("nextcur", "i4", ())
            f.createVariable("rmin", "f8", ())
            f.createVariable("zmin", "f8", ())
            f.createVariable("rmax", "f8", ())
            f.createVariable("zmax", "f8", ())

            f.createVariable(
                "coil_group", "S1", ("external_coil_groups", "stringsize")
            )
            f.createVariable("mgrid_mode", "S1", ("dim_00001",))
            f.createVariable("raw_coil_cur", "f8", ("external_coils",))

            # Add br_001, bp_001, bz_001, br_002, bp_002, bz_002, etc.
            # and keep array in Fortran ordering (ir, jz, kp) -> (kp, jz, ir)
            for kgroup, group in enumerate(self.coil_group):
                f.createVariable(
                    f"br_{kgroup+1:03}", "f8", ("phi", "zee", "rad")
                )
                f.createVariable(
                    f"bp_{kgroup+1:03}", "f8", ("phi", "zee", "rad")
                )
                f.createVariable(
                    f"bz_{kgroup+1:03}", "f8", ("phi", "zee", "rad")
                )

            # Add data
            f.variables["ir"][:] = self.ir
            f.variables["jz"][:] = self.jz
            f.variables["kp"][:] = self.kp
            f.variables["nfp"][:] = self.nfp
            f.variables["nextcur"][:] = self.coil_current.size
            f.variables["rmin"][:] = self.rmin
            f.variables["zmin"][:] = self.zmin
            f.variables["rmax"][:] = self.rmax
            f.variables["zmax"][:] = self.zmax

            f.variables["coil_group"][:] = self.coil_group
            f.variables["mgrid_mode"][:] = self.mgrid_mode
            f.variables["raw_coil_cur"][:] = self.coil_current

            for kgroup, group in enumerate(self.coil_group):
                f.variables[f"br_{kgroup+1:03}"][:] = self.br[kgroup]
                f.variables[f"bp_{kgroup+1:03}"][:] = self.bp[kgroup]
                f.variables[f"bz_{kgroup+1:03}"][:] = self.bz[kgroup]
