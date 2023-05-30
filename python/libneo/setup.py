from numpy.distutils.core import setup, Extension

setup(
    name='libneo',
    version='0.0.1',
    description='libneo',
    author='Christopher Albert',
    ext_modules=[
        Extension('magfie',
                  sources=['f2py_magfie.f90',
                           '../../src/magfie/magfie_vmec.f90',
                           '../../src/magfie/spline_vmec_data.f90'],
                  include_dirs=['../../build'],
                  libraries=['liblibneo.so'],
                  library_dirs=['../../build/'])
    ],
    data_files=[('libneo', ['../../build/liblibneo.so'])]
)
