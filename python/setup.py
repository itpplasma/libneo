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
