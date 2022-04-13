"""

    Get toroidal angle


"""

__author__='Matthias Willensdorfer (mwillens@ipp.mpg.de)'
__date__='21 Nov 2014'
__version__='1.0'

def DCN( channel = 'H_1' ):
    H_1 = 303.70 - 67.5
    H_2 = H_1
    H_3 = H_1
    H_4 = 304.70 - 67.5
    H_5 = H_1
    if (channel == 'H_1') | (channel == 'H_2') | (channel == 'H_3'):
        return H_1
    else:
        return H_5

def RPS( ):
    return 180.- 67.5

def LIB( ):
    return 260.64- 67.5

def ECE( ):
    return 191.25

def SXR( ):
    return 593.9 - 360.

def probes( ):
    return 225. - 67.5

def CTS( ):
    return 104.2

def VTN( ):
    return 123.90 - 67.5


def FRS( ):
    return 163.52 - 67.5

