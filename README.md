# libneo
Common code for plasma codes of ITPcp, in particular for different Neo-2 versions.

##aug
This is mostly python code for use in AUG intranet to extract data from the AUG shotfile database. The following libraries are used by this.

###eqtools_modified
Modified version of the [eqtools](https://github.com/PSFCPlasmaTools/eqtools) library that is available on GITHUB. Modifications where made because the code is not maintained anymore and it did not work in the beginning.

###wlsmodpy
Python code to write g-files from equilibria.

###mwilmodpy
Python code to include PSL response of coil files. There is way more in there but not used by now.

##doc

##matlab
Matlab library that contains various functions, classes and code interfaces for general and specific uses.

###BALANCE
Interface for the balance code. Contains Matlab code that is used to run the code, make pre- and post processing. Contains Fortran source files that are used in the interface.

###EFIT
A class to read/write/modify equilibrium g-files (commonly wrongly called efit files).

###GPEC_interface
A minimalist interface to run the code GPEC together with DCON.

###InputFile
A class that can be used to read/create/modify inputfiles that contain Fortran Namelists.

###KIN
A class that can write .kin files (kinetic profiles saved column wise)

###KiLCA_interface
Interface to prepare/run/postprocess the code KiLCA. A compiled version of KiLCA is needed (see [here](https://github.com/itpplasma/KiLCA))

###Shared
Shared classes and functions of the Matlab library.

###Utility
Classes and functions that have a general use for many application.

###mnDAT
Class to write mnDAT files.

##python

##src

##tests
