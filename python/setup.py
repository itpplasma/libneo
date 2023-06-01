import os
import setuptools
from numpy.distutils.core import setup, Extension

libneo_build_path = os.environ.get('LIBNEO_BUILD_PATH', '../build')

setup(
    name='libneo',
    version='0.0.1',
    description='libneo',
    author='Christopher Albert',
    packages=['libneo'],
    entry_points={
        'console_scripts': [
            'vmec_to_boozer = libneo.boozer:main',
            'convert_neo_in_to_fortran_namelist = '
                'libneo.convert_neo_in_to_fortran_namelist:main',
            'get_header_data_vmec = libneo.getHeaderDataVMEC:main'
        ],
    },
    ext_modules=[
        Extension('libneo.magfie',
                  sources=['f2py_interfaces/f2py_magfie.f90',
                           '../src/magfie/magfie_vmec.f90',
                           '../src/spline_vmec_data.f90'],
                  include_dirs=[libneo_build_path],
                  libraries=['libneo'],
                  library_dirs=[libneo_build_path])
    ],
)
