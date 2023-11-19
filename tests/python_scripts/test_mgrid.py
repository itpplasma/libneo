import numpy as np

from libneo.mgrid import MgridFile, list_to_bytearray, bytearray_to_list

def test_bytearray_to_list():
    bytearray = np.ma.masked_array(
        data=[[b'H', b'C', b'1', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
                [b'H', b'C', b'2', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
                [b'O', b'V', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
                [b'I', b'S', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
                [b'I', b'V', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ']],
        mask=False,
        fill_value=b'N/A',
        dtype='|S1')

    assert bytearray_to_list(bytearray) == ['HC1', 'HC2', 'OV', 'IS', 'IV']

def test_list_to_bytearray():
    input_strings = ['HC1', 'HC2', 'OV', 'IS', 'IV']
    expected = np.array(
        [[b'H', b'C', b'1', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
       [b'H', b'C', b'2', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
       [b'O', b'V', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
       [b'I', b'S', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' '],
       [b'I', b'V', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
        b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ']],
        dtype='|S1')

    assert np.array_equal(list_to_bytearray(input_strings), expected)


def test_mgrid_construct():
    f = MgridFile(
        ir=64, jz=64, kp=36, rmin=0.1, zmin=-1.0, rmax=1.0, zmax=1.0
    )
    assert f.ir == 64
    assert f.jz == 64
    assert f.kp == 36
    assert f.rmin == 0.1
    assert f.zmin == -1.0
    assert f.rmax == 1.0
    assert f.zmax == 1.0
    assert f.nfp == 1
    assert f.mgrid_mode == "S"
    assert f.coil_group == []
    assert f.coil_current == []
    assert f.br == []
    assert f.bp == []
    assert f.bz == []

def test_mgrid_read():
    f = MgridFile.from_file('/proj/plasma/DATA/LHD/albert/vmec/makegrid_alternative/mgrid_lhd_nfp10.nc')
    print(f)

def test_mgrid_read_write():
    f = MgridFile.from_file('/proj/plasma/DATA/LHD/albert/vmec/makegrid_alternative/mgrid_lhd_nfp10.nc')
    f.write('/tmp/mgrid_test.nc')


def test_mgrid_read_write_read():
    f = MgridFile.from_file('/proj/plasma/DATA/LHD/albert/vmec/makegrid_alternative/mgrid_lhd_nfp10.nc')
    f.write('/tmp/mgrid_test.nc')
    g = MgridFile.from_file('/tmp/mgrid_test.nc')

    # Assert if all variables are the same
    assert f.ir == g.ir
    assert f.jz == g.jz
    assert f.kp == g.kp
    assert f.nfp == g.nfp
    assert f.rmin == g.rmin
    assert f.zmin == g.zmin
    assert f.rmax == g.rmax
    assert f.zmax == g.zmax
    assert f.mgrid_mode == g.mgrid_mode
    assert np.array_equal(f.coil_group, g.coil_group)
    assert np.array_equal(f.coil_current, g.coil_current)
    assert np.array_equal(f.br, g.br)
    assert np.array_equal(f.bp, g.bp)
    assert np.array_equal(f.bz, g.bz)
