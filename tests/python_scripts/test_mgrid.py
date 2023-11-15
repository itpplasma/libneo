from libneo.mgrid import MgridFile

def test_mgrid_read():
    f = MgridFile('/proj/plasma/DATA/LHD/albert/vmec/makegrid_alternative/mgrid_lhd_nfp10.nc')
    print(f)


