"""
Handle mgrid data
"""
import textwrap
import numpy as np
from netCDF4 import Dataset


def list_to_bytearray(input_strings, max_length=30):
    data_list = [[char.encode() for char in string] for string in input_strings]

    # Pad the lists with empty bytes to make them of equal length
    data_list = [
        sublist + [b" "] * (max_length - len(sublist)) for sublist in data_list
    ]

    return np.array(data_list, dtype="S1")


def bytearray_to_list(input_bytearray):
    return [
        b"".join(bytearray).decode().strip() for bytearray in input_bytearray
    ]


class MgridFile:
    def __init__(
        self,
        ir,
        jz,
        kp,
        rmin,
        zmin,
        rmax,
        zmax,
        nfp=1,
        mgrid_mode="S",
        coil_group=None,
        coil_current=None,
        br=None,
        bp=None,
        bz=None,
    ):
        # Set self.ir, self.jz, etc. dynamically
        vars(self).update((k, v) for k, v in vars().items() if k != "self")

        if self.coil_group is None:
            self.coil_group = []
        if self.coil_current is None:
            self.coil_current = []
        if self.br is None:
            self.br = []
        if self.bp is None:
            self.bp = []
        if self.bz is None:
            self.bz = []

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

    def __repr__(self):
        return self.__str__()

    @classmethod
    def from_file(cls, filename):
        with Dataset(filename, "r") as f:
            coil_group = bytearray_to_list(f.variables["coil_group"][:])
            coil_current = f.variables["raw_coil_cur"][:]
            ir = f.dimensions["rad"].size
            jz = f.dimensions["zee"].size
            kp = f.dimensions["phi"].size

            nfp = f.variables["nfp"][:]
            rmin = f.variables["rmin"][:]
            zmin = f.variables["zmin"][:]
            rmax = f.variables["rmax"][:]
            zmax = f.variables["zmax"][:]
            mgrid_mode = f.variables["mgrid_mode"][:]

            br = []
            bp = []
            bz = []
            for kgroup, group in enumerate(coil_group):
                br.append(f.variables[f"br_{kgroup+1:03}"][:])
                bp.append(f.variables[f"bp_{kgroup+1:03}"][:])
                bz.append(f.variables[f"bz_{kgroup+1:03}"][:])

        return cls(
            ir=ir,
            jz=jz,
            kp=kp,
            rmin=rmin,
            zmin=zmin,
            rmax=rmax,
            zmax=zmax,
            nfp=nfp,
            mgrid_mode=mgrid_mode,
            coil_group=coil_group,
            coil_current=coil_current,
            br=br,
            bp=bp,
            bz=bz,
        )

    def write(self, filename):
        with Dataset(filename, "w") as f:
            # Add dimensions
            f.createDimension("stringsize", 30)
            f.createDimension("external_coil_groups", len(self.coil_group))
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

            f.variables["coil_group"][:] = list_to_bytearray(
                self.coil_group, 30
            )
            f.variables["mgrid_mode"][:] = self.mgrid_mode
            f.variables["raw_coil_cur"][:] = self.coil_current

            for kgroup, group in enumerate(self.coil_group):
                f.variables[f"br_{kgroup+1:03}"][:] = self.br[kgroup]
                f.variables[f"bp_{kgroup+1:03}"][:] = self.bp[kgroup]
                f.variables[f"bz_{kgroup+1:03}"][:] = self.bz[kgroup]

    def add_coil_group(self, group_name):
        self.coil_group.append(group_name)

    def add_coil(self, group_name, current):
        if not group_name in self.coil_group:
            self.coil_group.append(group_name)
        self.coil_current.append(current)
