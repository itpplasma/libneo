import setuptools
from numpy.distutils.core import setup, Extension

setup(
    name='libneo',
    version='0.0.1',
    description='libneo',
    author='Christopher Albert',
    packages=['libneo'],
    ext_modules=[
        Extension('libneo.magfie',
                  sources=['f2py_magfie.f90',
                           '../../src/magfie/magfie_vmec.f90',
                           '../../src/spline_vmec_data.f90'],
                  include_dirs=['../../build'],
                  libraries=['libneo'],
                  library_dirs=['../../build/'])
    ],
)
