# For similar setup, see https://github.com/simonsobs/pspy/blob/master/setup.py

import os
import shutil
import subprocess
import sys

from pathlib import Path
from setuptools import setup, find_packages
from numpy.distutils.command import build_ext
from numpy.distutils.core import Extension

name = "libneo"
libneo_build_path = os.environ.get("LIBNEO_BUILD_PATH", "../../../build")


class MakeBuild(build_ext):

    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        print("Building extension", ext.name)

        self.build_path = Path(self.build_temp)
        os.environ["PYTHON"] = sys.executable

        if self.build_path.is_dir(): shutil.rmtree(self.build_temp)
        shutil.copytree("../src/", self.build_temp)
        shutil.copy(".f2py_f2cmap", self.build_temp)

        sources = " ".join(ext.sources)
        options = self._options_if_defined(ext)
        cmd = f"f90wrap -k .f2py_f2cmap -m {ext.name} {sources} {options}"
        print(cmd)
        res = subprocess.run(cmd, cwd=self.build_temp, shell=True)
        if res.returncode != 0:
            raise Exception("f90wrap failed")

        stripped_sources = [s.split("/")[-1] for s in ext.sources]
        allfiles = os.listdir(self.build_temp)

        f90wrap_sources = [f for f in allfiles if f.startswith("f90wrap_") if f.split("f90wrap_")[1] in stripped_sources or f == "f90wrap_toplevel.f90"]

        other_f90_sources = [f for f in allfiles if f.endswith(".f90") if f in stripped_sources if "f90wrap_"+f not in f90wrap_sources]
        other_sources =  [f for f in stripped_sources if f.endswith(".f")]

        f90wrap_sources = " ".join(f90wrap_sources)
        other_f90_sources = " ".join(other_f90_sources)
        other_sources = " ".join(other_sources)

        include_dirs = " ".join([f"-I{d}" for d in ext.include_dirs])
        library_dirs = " ".join([f"-L{d}" for d in ext.library_dirs])
        libraries = " ".join([f"-l{d}" for d in ext.libraries])

        cmd = f"gfortran -fPIC -c {' '.join(ext.sources)} {include_dirs}"
        print(cmd)
        res = subprocess.run(cmd, cwd=self.build_temp, shell=True)
        if res.returncode != 0:
            raise Exception("Compilation failed")
        objects = " ".join([f.split(".f")[0]+".o" for f in stripped_sources])

        cmd = f"f2py-f90wrap -c -m _{ext.name} {objects} {f90wrap_sources} {include_dirs} {library_dirs} {libraries}"
        print(cmd)
        res = subprocess.run(cmd, cwd=self.build_temp, shell=True)
        if res.returncode != 0:
            raise Exception("f2py failed")

        if not os.path.exists(self.build_lib): os.makedirs(self.build_lib)
        print(self.build_lib)


        outdir_lib = os.path.join(self.build_lib, "_" + name)
        outdir_python = os.path.join(self.build_lib, name + "/")

        if not os.path.exists(outdir_lib): os.makedirs(outdir_lib)
        if not os.path.exists(outdir_python): os.makedirs(outdir_python)

        for f in self.build_path.glob('*.so'):
            shutil.copy2(f, outdir_lib)
        for f in self.build_path.glob('*.py'):
            shutil.copy2(f, outdir_python)





    def _options_if_defined(self, ext):
        print("EXT: ", ext.__dict__)
        if not "f2py_options" in ext.__dict__:
            return ""
        skip = self._skip_if_defined(ext.f2py_options)
        return skip

    def _skip_if_defined(self, options):
        if not "skip:" in options:
            return ""
        index = options.index("skip:")
        next_colon = options[index+1:].index(":")
        skip_options = options[index+1:index+1+next_colon]
        return "--skip " + " ".join(skip_options)

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
        Extension('magfie',
                  sources=['f2py_interfaces/f2py_magfie.f90',
                           'magfie/magfie_vmec.f90',
                           'spline_vmec_data.f90'],
                  include_dirs=[libneo_build_path],
                  libraries=['libneo'],
                  library_dirs=[libneo_build_path]),
        Extension('efit_to_boozer',
                sources=[
                     'f2py_interfaces/f2py_efit_to_boozer.f90',
                     'efit_to_boozer/SRC/odeint_allroutines.f',
                     'efit_to_boozer/SRC/efit_to_boozer.f90',
                     'efit_to_boozer/SRC/field_divB0.f90',
                     'efit_to_boozer/SRC/efit_to_boozer_mod.f90',
                     'efit_to_boozer/SRC/spline5_RZ.f90',
                     'efit_to_boozer/SRC/spl_three_to_five_mod.f90',
                     'efit_to_boozer/SRC/bdivfree_coul.f90',
                     'efit_to_boozer/SRC/field_line_integration_for_Boozer.f90',
                     'efit_to_boozer/SRC/plag_coeff.f90',
                     'efit_to_boozer/SRC/binsrc.f90',
                     'efit_to_boozer/SRC/rhs.f90',
                     'efit_to_boozer/SRC/spline_and_interpolate_magdata.f90'
                     ],
                    f2py_options = [
                        'skip:','oddorderspline', 'polleg',
                        'binomial', 'odeint_allroutines', 'alloc_odeint',
                        'odeint', 'rkck', 'rkqs', 'binsrc', 'spl_five_reg', ':'
                    ],
                    include_dirs=['efit_to_boozer/OBJS'],
                    libraries=['blas', 'lapack'],
                ),
        Extension('interpolate',
                    sources=[
                        'interpolate.f90',
                    ],
                ),
    ],
    cmdclass = {"build_ext" : MakeBuild},
)
