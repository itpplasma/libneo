"""
Handle mgrid data
"""
from netCDF4 import Dataset

class MgridFile:
    def __init__(self, filename=None):
        if filename:
            self.read(filename)

    def __str__(self):
        return "mgrid"
    
    def read(self, filename):
        with Dataset(filename, "r") as f:
            pass


    def write(self, filename):
        f = Dataset(filename, "w")

        # Add dimensions
        f.createDimension("stringsize", 30)
        f.createDimension("external_coil_groups", self.coil_group.size)
        f.createDimension("dim_00001", 1)
        f.createDimension("external_coils", self.coils.size)
        f.createDimension("rad", self.rad.size)
        f.createDimension("zee", self.zee.size)
        f.createDimension("phi", self.phi.size)

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

        f.createVariable("coil_group", "S1", ("external_coil_groups","stringsize"))
        f.createVariable("mgrid_mode", "S1", ("dim_00001",))
        f.createVariable("raw_coil_cur", "f8", ("external_coils",))

        # TODO: Add br_001, bp_001, bz_001, br_002, bp_002, bz_002, etc.
