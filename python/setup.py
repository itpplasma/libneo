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
    install_requires=[
        'numpy'
    ],
    ext_modules=[
        Extension('libneo.magfie',
                  sources=['f2py_interfaces/f2py_magfie.f90',
                           '../src/magfie/magfie_vmec.f90',
                           '../src/spline_vmec_data.f90'],
                  include_dirs=[libneo_build_path],
                  libraries=['libneo'],
                  library_dirs=[libneo_build_path]),
        Extension('libneo.efit_to_boozer',
                sources=[
                     '../../efit_to_boozer/SRC/efit_to_boozer.f90',
                     '../../efit_to_boozer/SRC/field_divB0.f90',
                     '../../efit_to_boozer/SRC/efit_to_boozer_mod.f90',
                     '../../efit_to_boozer/SRC/spline5_RZ.f90',
                     '../../efit_to_boozer/SRC/spl_three_to_five_mod.f90',
                     '../../efit_to_boozer/SRC/bdivfree_coul.f90',
                     '../../efit_to_boozer/SRC/field_line_integration_for_Boozer.f90',
                     '../../efit_to_boozer/SRC/plag_coeff.f90',
                     '../../efit_to_boozer/SRC/binsrc.f90',
                     '../../efit_to_boozer/SRC/rhs.f90',
                     '../../efit_to_boozer/SRC/spline_and_interpolate_magdata.f90',
                     'f2py_interfaces/f2py_efit_to_boozer.f90'
                     ],
                  include_dirs=[libneo_build_path],
                  libraries=['libneo'],
                  library_dirs=[libneo_build_path])
    ]
)
