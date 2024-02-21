# File formats of Magnetic field data

In this section we explain the different file formats that one might encounter when dealing with magnetic field data. We state in which context they usually appear, as well as codes generating it/converting it to other formats. 

For syntax: The name of each section is given by the file format (i.e. how the data is structured/arranged in the file). In paranthese one finds the type of field (2D or 3D) that is usually stored, as well as a common code that generates the files in question. Finally the file type (different from format) i.e. the `.extension` of the file (e.g. textfile, binary netcdf, ...) is listed as well.

## gfiles/EQDSK format (axisymmetric 2D equillibrium, EFIT, textfile)

[Original documentation from 1997](https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf) (note that `xdum` stands for a place holder to keep the formating of columns) <br>

Units: SI (**note**: converters using gfiles as source may give out the quantities in cgs units after conversion) <br>
entries per line: 5 (or less if end of points) <br>
values format: ~ 16 characters (including sign and exponent notation) <br>

While the general structure (see link above) should be the same for all files of that format, additional linebreaks and different formating of floating point numbers (trailing spaces) makes reading these textfiles with practically most programming languages besides FORTRAN cumbersome. Therefore, if one works e.g. in Python it not possible to have one uniform reader routine for these textfiles. The user still has to choose the right reader routine for their present case, modify them accordingly or write their own. For FORTRAN however, this should not be a problem. If however, FORTRAN should actually fail to read it, a preprocessing of the file should be prefered, instead of including an extra reader for these edgecases in the main FORTRAN code.

To give some orientation, for files which look like

     Equil_2021_PMI_QH_mode_betap_1d04_li_1d02_Ip_18 3 360 400
     18 20 8.9344 0 0
     9.467 0.0583 302.6816 135.5856 -5.7463                1st block of values
     -18270000.0013 -167.096 0 9.467 0       <- shot info           |
     0.0583 -167.096 0 9.467 0                                      V
     -5.09312766e+01 -5.09337203e+01 -5.09361562e+01 -5.09385842e+01 -5.09410043e+01 
     -5.09434164e+01 -5.09458207e+01 -5.09482171e+01 -5.09506056e+01 -5.09529863e+01 
    ...
     -5.13398373e+01 -5.13398485e+01 -5.13398619e+01 -5.13398682e+01 -5.13398665e+01 
                                                                                    <-\n
     +1.92499100e+06 +1.91352411e+06 +1.90209288e+06 +1.89069842e+06 +1.87934049e+06 
    ...

where there is a line break in between the data blocks and the shot information (width in $R$ and $Z$ of box where data points are contained, leftmost $R$, ...) are listed in a different format than the quantity values in the subsequent blocks, the [eqdsk.py](https://github.com/itpplasma/libneo/blob/main/python/libneo/eqdsk.py) Python class from libneo, can usually be used. You have to add libneo/python/ to your library [search path](https://note.nkmk.me/en/python-import-module-search-path/), see also [Libraries](https://gitlab.tugraz.at/plasma/info/-/tree/main/TOPIC_Compilation_and_Libraries/Libraries.md). Alternatively, you can just make a symbolic link in the same directory you want to import the reader.

    from libneo import eqdsk # You have to add libneo/python/ to your library search path
                             # or have symbolic link in same directory
    eqdsk_data = eqdsk.eqdsk_file(filename) # Creates a python dictionary eqdsk_data

Note that this reader routine was designed with having additional coil information positioned after the contents described in the linked documentation

A reader that is more simpler and oriented on files that look like

          FROM CHEASE BUT COCOS=02 ,SI UNITS20220422   0 181 181                     ---
     0.633506867E+01 0.125330640E+02 0.893800000E+01 0.577127254E+01 0.866000590E-01  ^
     0.948609356E+01 0.866000590E-01-0.211906390E+02 0.000000000E+00 0.574000000E+01 shot
     0.180296699E+08-0.211906390E+02 0.000000000E+00 0.948609356E+01 0.000000000E+00 info
     0.866000590E-01 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 _V_
     0.505304811E+02 0.505557945E+02 0.505912263E+02 0.506328321E+02 0.506773378E+02  ^
     0.507225815E+02 0.507671753E+02 0.508101542E+02 0.508509820E+02 0.508893147E+02  |
     0.509249999E+02 0.509580936E+02 0.509885957E+02 0.510166927E+02 0.510425551E+02 1st
    ...                                                                             block
     0.512468684E+02 0.512563257E+02 0.512673536E+02 0.512793829E+02 0.512918161E+02  |
     0.513041200E+02                                                                 _V_
     0.153492648E+07 0.151262906E+07 0.148577727E+07 0.145611416E+07 0.142509588E+07
    ...

where there is no additional linebreaks and the shot data is formated in the same ways as the actual data, readers like [poor_eqdsk.py](https://github.com/itpplasma/ntv-demo/blob/main/python/poor_eqdsk.py) from the ntv-demo repo may be used. However, again it has to be stressed that there is not one uniform reader and you might very well need to write one yourself. We will therefore explain the parts of the file, so that you can extract the necessary information yourself if needed.

- `PsiaxisVs` ... the **poloidal** magnetic flux at the magnetic axis **normalized** by a factor $1/(2\pi)$ <br>
to get the actual poloidal flux (e.g. for the input of `VMEC`) one has to multiply by $2\pi$ first <br>
IF this value is $\neq 0$ it is the poloidal disk flux $\psi_\text{pol}^\text{disk}$, see also [lookup for flux coordinates](https://gitlab.tugraz.at/plasma/info/-/tree/main/TOPIC_MHD_Equilibria/Lookup_flux_coordinates.md) or D'haeseleer.
- `PsiedgeVs` ... the **poloidal** magnetic flux at the edge i.e. boundary flux surface **normalized** by a factor $1/(2\pi)$ <br>
to get the actual poloidal flux (e.g. for the input of `VMEC`) one has to multiply by $2\pi$ first

**Note** that depending on the further application you have to convert this into the **toroidal** flux (e.g. for `VMEC`) instead. This can be done with by using the `qprof` safety factor profil (see also the [FluxConverter](https://github.com/itpplasma/libneo/blob/main/python/libneo/flux_converter.py) of libneo) 

### TODO

The last two numbers specified in the header (e.g. in example above the two `181`) give the number of equidistant grid points in $R$- and $Z$-direction inside the box specified in the shot info below. 

$N_R$ ...  number of equidistant points in $R$-direction, second to last number in header <br>
$Z_R$ ...  number of points in $Z$-direction, last number in header <br>

There is very important point to keep in mind. The second to last number does not only specify $N_R$, but also

!!! <br>
$N_s$ ... number of equidistant points in $s$-direction, second to last number in header <br>
!!! <br>

It turns out that the majority of quantities provided in the EQDSK format (pressure, rotational transform, ...) are given not on the above specified $(R,Z)$ rectangular grid, but on an uniform poloidal flux grid (i.e. on the flux surfaces labeled by $s_\text{pol}$). Notable exception is the poloidal flux itself, which is specified on the rectangular grid.

![Visualization of some content of the eqdsk file](/doc/ExtraDocuments/pictures/eqdsk_data_example_1.png)

#### Boundary

defined by in terms of geometrical poloidal angle $\theta$ equidistant points, forming a closed flux surface

## MARSQ_file (axisymmetric 2D equillibrium, MARSQ, matlab-strcture)
load by MATLAB internal routine `load()` or for example Python `loadmat()` <br>
structure with field names like `R_MARS` or `rho_Boozer` (different coordinate systems possible/mixed) <br>
e.g. <br>

$\rho_{Boozer} = [1 \times 100]$ - equidistant array of BOOZER flux label <br>
$\vartheta_{Boozer} = [1 \times 50] $ - equidistant arra of BOOZER theta angle<br>
$R_{BOOZE} = [100 \times 50] $ &rarr; points in CYLINDRICAL coordinates ($R$) on a BOOZER equidistant grid<br>
$Z_{BOOZE} = [100 \times 50] $ &rarr; points in CYLINDRICAL coordinates ($Z$) on a BOOZER equidistant grid<br>
$Bz_{Boozer} = [100 \times 50] $ &rarr; CYLINDRICAL covariant component on a BOOZER equidistant grid <br>

flux label is actually stored as $\rho_{Boozer} =\sqrt{s}$ i.e. square root of normalised, poloidal flux.

TODO

## XPLASMA.OUT (,MARS,textfile)
TODO

## vmec_files (3D equillibrium, `VMEC`, binary netcdf file)
TODO

## boozX_files (3D equillibrium, `BOOZX`, binary netcdf file)
TODO

## `.bc`-files or boozerfiles (2D or 3D equi., ITP converter routines, textfile)
TODO

The length of the values was shortend to the usual precision for better readability.

    CC Boozer-coordinate data file
    CC Version: 
    CC Author: ert
    CC shot:    0
    m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]
    24     0  1000    1  -3.239690e+02   4.075530e+00   9.486130e+00
    s            iota        Jpol/nper       Itor         pprime      sqrt (0,0)
                                [A]           [A]          [Pa]      (dV/ds)/nper
    5.0000e-04   9.1687e-01   2.5652e+08   2.6486e+04  -4.4077e+06  -3.5691e+03
    m   n   rmnc [m]  rmns [m]  zmnc [m]  zmns [m]  vmnc [ ]  vmns [ ]  bmnc [T]  bmns [T]
    0   0   9.4e+00   0.0e+00   8.6e-02   0.0e+00  -5.6e-08   0.0e+00   5.4e+00   0.0e+00
    1   0   8.1e-02  -9.3e-04  -5.8e-05   1.1e-01  -1.0e-08  -4.8e-07  -4.6e-02   5.3e-04
    2   0   1.1e-03  -2.4e-05   2.2e-05   1.3e-03   6.6e-08  -2.7e-06  -3.6e-04   1.1e-05
    ...

### Meaning of quantities

Here we list the meaning of the quantities in the file. While the actual names may differ between versions of Boozer files, their position in the file should stay the same, which is the relevant criterion for readers implemented in our codes like [NEO-2](https://gitlab.tugraz.at/plasma/info/-/tree/main/codes/NEO2.md).

#### Quantities common for all surfaces
- `m0b` & `n0b` <br>
Give (or should give) the maximum mode number in poloidal ($m$) or toroidal ($n$) regard. Therefore in the file above there are $24 + 1$ different poloidal mode numbers and $0 + 1$ toroidal mode numbers. Depending on the specific form of the Boozer file, one can calculate the total number of modes present per surface (see below).

- `nsurf` <br>
Number of surfaces in the file

- `nper` <br>
Number of period, also often called `nfp`, see [Number of Field Periods](/doc/ExtraDocuments/Conventions.md#number-of-field-periods) in the Convention section.

- `flux` <br>
The toroidal magnetic flux at the edge/seperatrix (NOT normalised by a factor $2\pi$). The flux has a sign, therefore depending on the coordinate system used in the Boozer file & in the application reading the Boozer file, an additional conversion has to be done.

- `a` <br> TODO
- `R` <br> TODO

#### Quantities specific to each surface

TODO
### Different formats

While the file presented above is a common variant at our institute, it is by far not the only version. Besides different names for the quantities, there are other differences like number of columns or quantities missing and replaced by others.

#### The form of the Fourier series

The Fourier series is expressed in terms of real coefficients i.e. the series consists of $\sin$ and cosine instead of complex exponentials. This is signified by the ending of the coefficients i.e.

- `*s` for coefficients of the $\sin$ part (e.g. `zmns`) and
- `*c` for coefficients of the cosine part (e.g. `rmnc`).

A very general way of writing such a Fourier series is

$\sum_{m=-\infty}^{\infty}\sum_{n=-\infty}^{\infty}f^\text{c}_{mn}\cos{(m\theta + n\varphi)} + f^\text{s}_{mn}\sin{(m\theta + n\varphi)}$.

The respective coefficients are connected to the complex coefficients of the complex series

$\sum_{m=-\infty}^{\infty}\sum_{n=-\infty}^{\infty}f_{mn}\exp{(i(m\theta + n\varphi))}$

by

$f^\text{s}_{mn} = -\Im{f_{mn}}$ <br>
$f^\text{c}_{mn} = \Re{f_{mn}}$.

As now for real functions like physical quantities, the complex coefficients have to fulfill

${f_{mn}}^* = f_{-m-n}$,

which further leads to

$f^\text{s}_{mn} = -f^\text{s}_{-m-n}$ <br>
$f^\text{c}_{mn} = +f^\text{c}_{-m-n}$

One therfore only needs to know half of the coefficients. This means, that one can choose **either** $m$ or $n$, to be non-negative (you only can do this for one of them, as else you would be missing modes). We commonly choose the poloidal mode $m$ for that. The coefficients are then redefined as

${f^\text{s}_{00}}' = 0$<br>
$f^\text{s}_{mn}\sin{(m\theta + n\varphi)} + f^\text{s}_{-m-n}\sin{(-m\theta - n\varphi)}=2f^\text{s}_{mn}\sin{(m\theta + n\varphi)}$ <br>
&rarr;<br>
${f^\text{s}_{mn}}' \equiv 2f^\text{s}_{mn}$<br>

${f^\text{c}_{00}}' = f^\text{c}_{00}$<br>
$f^\text{c}_{mn}\cos{(m\theta + n\varphi)} + f^\text{c}_{-m-n}\cos{(-m\theta - n\varphi)}=2f^\text{c}_{mn}\cos{(m\theta + n\varphi)}$ <br>
&rarr;<br>
${f^\text{c}_{mn}}' \equiv 2f^\text{c}_{mn}$<br>

and the sum goes like

$\sum_{m=0}^{\infty}\sum_{n=-\infty}^{\infty}{f^\text{c}_{mn}}'\cos{(m\theta + n\varphi)} + {f^\text{s}_{mn}}'\sin{(m\theta + n\varphi)}$

Alternatively, one can go one step further. As

${f^\text{s}_{0n}}' = 2{f^\text{s}_{0n}} = -2f^\text{s}_{-0-n} = -2f^\text{s}_{0-n} = -{f^\text{s}_{0-n}}'$ <br>
${f^\text{c}_{0n}}' = 2{f^\text{c}_{0n}} = +2f^\text{c}_{-0-n} = +2f^\text{c}_{0-n} = +{f^\text{c}_{0-n}}'$

one can again recombine to (analoge to the previous one)

${f^\text{s}_{0n}}'' \equiv 2{f^\text{s}_{0n}}'$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for $n>0$<br>
${f^\text{s}_{0n}}'' \equiv 0$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for $n<0$ &rarr; left out in Boozer file<br>
${f^\text{s}_{00}}'' \equiv {f^\text{s}_{00}}' = 0$<br>
${f^\text{s}_{mn}}'' \equiv {f^\text{s}_{mn}}'$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for $m \ne 0$

${f^\text{c}_{0n}}'' \equiv 2{f^\text{c}_{0n}}'$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for $n>0$<br>
${f^\text{c}_{0n}}'' \equiv 0$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for $n<0$ &rarr; left out in Boozer file<br>
${f^\text{c}_{00}}'' \equiv {f^\text{c}_{00}}' = {f^\text{c}_{00}}$<br>
${f^\text{c}_{mn}}'' \equiv {f^\text{c}_{mn}}'$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for $m \ne 0$

In this new definition, both for $\sin$ and cosine ,there are no negative modes if $m=0$. One can now choose (but does not have to) leave out these zeros when writing the Boozer file. Depening on the choice, the number of present modes in the file changes (see below). 

However, for the reconstruction the definition of the modes is actually not relevant (assuming one can read out all the rows properly from the file). Because no matter which definition of coefficient is present, as long as one adds all **present** modes up with the correct mode numbers $m$ & $n$, one gets the correct quantity back.

TD;LR: While from the understanding there are different definitions of the fouriercoefficients, as long as one reads the **present** modes out literally (i.e.

`m=0 n=-1 rmnc=0.5 rmns=0.2` &rarr; $0.5\cos{(0\theta - \varphi)} + 0.2\sin{(0\theta -\varphi)}$<br>
`m=1 n=+1 rmnc=0.2 rmns=0.5` &rarr; $0.2\cos{(\theta + \varphi)} + 0.5\sin{(\theta +\varphi)}$<br>
`m=0 n=-1 rmnc=0.0 rmns=0.0` &rarr; $0.0\cos{(0\theta - \varphi)} + 0.0\sin{(0\theta -\varphi)}$<br>
$\ldots$

) and adds them all up, one does always end up with the correct value.

#### number of modes per surface

using the above `m0b` and `n0b` one can deduce how many modes per surface there are, when considering how the modes are listed. Common ways are

- $m \ge 0$ and for each $m \neq 0$ toroidal modes $n \in [-\text{n0b},\text{n0b}]$, for $m=0 \rightarrow n \in [0,\text{n0b}]$ <br>
The total number of modes is therefore $(\text{n0b}+1) + \text{m0b}*(2*\text{n0b}+1)$ <br>
This format is sometimes refered to at ITP as "Strumberger Format"

- $m \ge 0$ and for each $m$ toroidal modes $n \in [-\text{n0b},\text{n0b}]$ <br>
The total number of modes is therefore $(\text{m0b}+1)*(2*\text{n0b}+1)$

#### Assuming stellarator symmetry

TODO

If stellarator symmetry is present, only the $\sin$ OR cosine contribution for a quantity e.g. the Radius $R$ is non-zero. Each quantity therefore then only has one column of coefficients assigned to it, reducing the overall columns to 6 instead of 8.