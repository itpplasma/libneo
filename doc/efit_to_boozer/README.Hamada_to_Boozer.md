# MARS Hamada to Boozer

Example run is in the directory `/proj/plasma/DATA/DEMO/BOOZER/example_conversion_mars_hamada_to_boozer`

First, run in this directory `efit_to_boozer.x` converter. Necessary inputs for this are

* `efit_to_boozer.inp`
* `field_divB0.inp`
* `convexwall.dat`
* an EQDSK file like `EQDSK_DEMO2019_q1_COCOS_02.OUT` in the above example folder

Last file here is EFIT file mentioned in `field_divB0.inp` and `convexwall.dat` should match this input file.

Second, run in the same directory `hamada_to_boozer.x` converter. It uses the iput for `efit_to_boozer.x`,
output from efit_to_boozer.x, and the following extra input specific for `hamada_to_boozer.x`:

* `mars_data_dimensions_and_extra.in`
* `DATATORQNTV.OUT`
* `RMZM_F.OUT`
* `TORQUENTV.OUT`

Last three files are output files of MARS (names are hard-coded in the converter). First file contains:
* mpolmax    - maximum by module poloidal mode index for Hamada spectrum from MARS
* nrhopol    - number of radial grid points for MARS data
* theta_H0   - value of Hamada poloidal angle at the point where symmetry flux angle = 0 (set to zero normally)
* theta_B0   - value of Boozer poloidal angle at the point where symmetry flux angle = 0 (set to zero normally)
* ntor       - toroidal mode number

For these inputs, `hamada_to_boozer.x` produces two files:

* `pert_fake.bc`
* `pert.bc`

First file is fake Boozer file where Fourier amplitudes in Hamada angles are interpreted as ampitudes in Boozer angles
Second file is actual Boozer file where Fourier amplitudes in Boozer angles are computed from Fourier amplitudes in Hamada angles
Description of the conversion procedure is in Section 2 of `EFIT_to_Boozer_and_Hamada_converter.pdf`.
Currently code interprets Hamada spectrum as the one with summation over n>0 and real part taken (no factor two appears for Boozer
spectrum, in contrast to what is in the PDF file).

May be of use: there is Boozer file reader in source directory:

`SRC/read_boozer_file_for_plotting.f90`

Complite it with gfortran (it needs no other files), link the perturbation Boozer file (`pert.bc` or `pert_fake.bc`) to

`booz_test.bc`

and run the executable. It reads coise and sine Fourier amlitudes from Boozer file and writes them to two files,
`bmnc.dat` and `bmns.dat` where the first column is s-value and the rest columns are 

$B_{mn}^c(mode)$ or $B_{mn}^s(mode)$ 

where "mode" is general Fourier mode index. First line in the output files starts with # (coment for gnuplot) and contains
the poloidal mode numbers as functions of general Fourier mode inex, m(mode). For Matlab plotting delete this line (or
comment it with %).
