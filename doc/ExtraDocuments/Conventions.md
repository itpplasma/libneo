# Conventions

This page will deal with explaining/linking the different sign conventions found in the `VMEC`, `BOOZX` and the ITP Plasma internal code base (i.e. the homemade converters).

## The gist of it

### Coordinate Systems

#### Cylindrical Coordinates
- radial coordinate $R$
- toroidal geometric angle $\varphi$
- height/symmetry axis coordinates $Z$
- most commonly right handed $(R,\varphi,Z)$, but left handed $(R,Z,\varphi)$ are also used

#### `VMEC` (custom system, non straight field lines)
- flux label/radius: $\Psi_V$ (= toroidal flux $\Psi_{\text{tor}}$)
- poloidal angle: $\theta_V$ (=/= geometric $\theta$)
- toroidal angle: $\varphi_V$ (= geometric $\varphi$)
- **Note** that because of the periodicity over the field periods, the actual angle used <br>
for the toroidal direction internally in `VMEC` is $\zeta_V \equiv N_\text{FP}*\varphi_V$ (see section about field periods below)
- left handed system $(\Psi_V,\theta_V,\varphi_V)$

#### `BOOZX` (Boozer coordinates, straight field lines, constant $B_{\theta_B}$ & $B_{\varphi_B}$ on a flux surface)
- flux label/radius: $\Psi_B = \Psi_V$
- poloidal angle: $\theta_B = \theta_V + \lambda(\Psi_V,\theta_V, \varphi_V) + G(\Psi_V,\theta_V, \varphi_V)$
- toroidal angle: $\varphi_B = \varphi_V + qG(\Psi_V,\theta_V, \varphi_V)$
- analogous to `VMEC` above an corresponding $\zeta_B \equiv N_\text{FP}*\varphi_V$ is used internally
- $\lambda$ is a transformation function to the straight fieldline version of  `VMEC` coordinates (see below)
- $p$ is a periodic function in the angles that transforms the straight fieldline coordinates into actual Boozer coordinates and the inverse rotational transform $q=1/\iota$ (intrinsic, physical quantity of the magnetic field, independant of coordinate system)
- left handed system $(\Psi_B,\theta_B,\varphi_B)$

#### straight fieldline `VMEC` (mentioned for completness, but NOT necessary/needed for workflow)
- flux label/radius: ${\Psi_V}^* = \Psi_V$ ($=\Psi_{\text{tor}}$) 
- poloidal angle: ${\theta_V}^* = \theta_V + \lambda(\Psi_V,\theta_V, \varphi_V) $
- toroidal angle: ${\varphi_V}^* = \varphi_V$ (= geometric $\varphi$)
- left handed system $({\Psi_V}^*,{\theta_V}^*,{\varphi_V}^*)$ 
- NOT used by any of the codes, but necessary if need straight field line system in postprocessing of `VMEC` data (without using `BOOZX` or ITP converters)

#### ITP Boozer coordinates = `BOOZX` coordinates
- BUT with convention to use the normalised flux as flux label/radius <br>
$\Psi_{\text{ITP}} = s = \Psi_{\text{tor}}/\Psi_{\text{tor}}(\text{Seperatrix}) = \Psi_{B}/\Psi_{B}(\text{Seperatrix})$
- **Attention**:The matter of $\zeta$ is different compared to the above cases (see Field period section below) 
- **Note** that using the normalised flux in the other systems can be also done without changing the distinct properties (therefore in the following, we use generically $s$ when talking about the radial coordinate)
- **Warning**: the choice of left/right handed system is not confirmed at this point TODO

### Sign of Kernel

When representing a quantity in terms of a spectral decomposition of two angle variables $(u,v)$ as

$f(s,u,v) = \sum_m\sum_n f_{mn}(s)\exp{(i*K(u,v,m,n))}$

the different codes apply a different definition of the Kernel $K = K(u,v,m,n)$:

| System   | Kernel K                 | Info |
|----------|--------------------------|------|
|  `VMEC`    |$(m\theta_V - n\varphi_V)$|see [`VMEC`-Wiki/Toroidal-Coordinate-Systems](https://princetonuniversity.github.io/STELLOPT/Toroidal%20Coordinates)
|  `BOOZX`   |$(m\theta_B - n\varphi_B)$|
|ITP Boozer|$(m\theta_B + n\varphi_B)$|

The convention of the Kernel has to be considered when using the modes $f_{mn}$ to reconstruct the quantity $f$ at a given point in the coordinate system.

<a id="realcoefficients"></a>
### Splitting of complex coefficients into sine and cosine contribution

In `VMEC` and such, the Fourier series is often written in terms of sine & cosine functions instead of the complex exponentials. 

$f(s,u,v) = \sum_m\sum_n f^c_{mn}(s)\cos{(K(u,v,m,n))} \pm f^s_{mn}(s)\sin{(K(u,v,m,n))}$

Notice the $\pm$ above. The complex coefficients $f_{mn}$ are split between these new, real coefficients $f^c_{mn}$ & $f^s_{mn}$, and depending on the choice of $\pm$ they have to be defined with an additional sign:

| System   | $f^c_{mn}$  | $f^s_{mn}$   | Info                              |
|----------|-------------|--------------|-----------------------------------|
|`VMEC`      |$\Re{f_{mn}}$|$-\Im{f_{mn}}$|using therefore $+$ for $\pm$ above|
|`BOOZX`     |$\Re{f_{mn}}$|$-\Im{f_{mn}}$|using therefore $+$ for $\pm$ above|
|ITP Boozer|$\Re{f_{mn}}$|$-\Im{f_{mn}}$|using therefore $+$ for $\pm$ above|

#### Further splitting (`VMEC`)

A praxis performed internally e.g. in `VMEC` (for example for the output of the magnetic axis), is a further splitting of harmonic functions into poloidal/toroidal parts. Using the trigonometric formulas one gets (using the `VMEC` Kernel)

$f(u,v) = \sum_{mn} f^\text{cc}_{mn}\cos{(m\theta_V)}\cos{(-n\varphi_V)} + f^\text{cs}_{mn}\cos{(m\theta_V)}\sin{(-n\varphi_V)} + f^\text{sc}_{mn}\sin{(m\theta_V)}\cos{(-n\varphi_V)} + f^\text{ss}_{mn}\sin{(m\theta_V)}\sin{(-n\varphi_V)}$.

We dropped the $s$ dependency here for readibility. Comparing with the above, one gets the coefficients straightforward as

$f^\text{cc}_{mn} = f^\text{ss}_{mn} = f^\text{c}_{mn}$ <br>
$f^\text{cs}_{mn} = f^\text{sc}_{mn} = f^\text{s}_{mn}$

One must keep in mind the form of the above series with this definition of the new coefficients. From the authors point of view it is rather weird that the negative sign is left inside the sine functions instead of defining the above coefficients with appropriated minus signs instead. However, despite that, this is the form that is used in `VMEC` internal representation (and assumed for the output) of certain quantities. An important example is the magnetic axis given back from `VMEC` in form of the coefficients [raxis_cc, z_axis_cs and such](https://gitlab.tugraz.at/plasma/info/-/blob/main/codes/`VMEC`.md#magneticaxis)

### Number of field periods

While a stellarator configuration is not symmetric in toroidal direction, it is however periodic. This means that after going a certain amount along the toroidal direction, the configuration repeats. Because of that the toroidal modes $n$ have to fulfill an additional periodicity criterion.

$f(u,v + t*2\pi/N_\text{FP}) = f(u,v)$ ... periodicity with period $2\pi/N_\text{FP}$ in toroidal direction

The condition implies that there are only toroidal mode numbers which are mutliples of $N_\text{FP}$, as one can check by applying the condition to an explicit series. This is achieved internally in the codes by using not an toroidal angle like $\varphi_V \in [0,2\pi[$ but rather $\zeta_V \equiv N_\text{FP}\varphi_V$ for the calculation. In the ouput of the codes this is however differently reflected

- `VMEC`: $n$ are saved as their actual form (so one has only mode-number $n$ that are multiples of `nfp`)
- `BOOZX`: same as `VMEC`
- ITP Boozer: the $n$ are stored as $n_\text{actual}/N_\text{FP}$, where $N_\text{FP}$ is called `nper` in the ITP Boozer files. They appear therefore not as multiples of the field period number but just like "-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,...". These saved modes have to be multipled by `nper` therfore again, before evaluating the modes in a Fourier series.

Note that for tokamaks, where there is no such concept of field periods, $N_\text{FP}$ is set to $1$.

## Left/right handed systems - The COCOS convention

[COCOS convention](https://crppwww.epfl.ch/~sauter/cocos/) tells you if you know the cocos number, the present orientation of coordinate systems (mostly tokamaks, axisymmetric codes, ...)

[Master thesis of Robert Babin](/TOPIC_MHD_Equilibria/sources/MA_Babin.pdf): see a good explanation of the nature of VMECs coordinate system (section 1.4/1.5)

As seen in the coordinate system section above, besides the choice of individual coordinates, the kind of system (left/right handed) can be different as well. One has to know the present kind, as the representation (the numbers we write down, have saved in e.g. output files) of the real world quantities (like a flux) are different. In short, one should be able from knowing the system type to imagine for example in which direction the poloidal flux flows and how to encode this into the different representations. In most cases it is merely a change of sign.

![Left & Right handed coordinate systems and the effect on represetation](/doc/ExtraDocuments/pictures/right_left_coordinate_systems.png)

TD;LR: You have to know which convention you are in, else you can not really do anything with the numbers you are reading out of the files.

### Example COCOS=2, interpreting EQDSK file content

Here is an example of an EQDSK file where the convention used is COCOS=2 (for this specific file). Using the associated coordinate system, one can reimagine the actual physical situation based on the "numbers" provided for the fluxes and safety factor. See also the quick [lookup](/doc/ExtraDocuments/Lookup_flux_coordinates.md) regarding straight field line coordinates and the safety factor (more in depth in D'heaseleer).

![COCOS2 example](/doc/ExtraDocuments/pictures/cocos2_physical_picture.png)


## Switching between Conventions

### Coordinate transformations to geometrical angle

The so called Boozer coordinates

$\theta_B = \theta_V + \lambda(\Psi_V,\theta_V, \varphi_V) + G(\Psi_V,\theta_V, \varphi_V)$ <br>
$\varphi_B = \varphi_V + qG(\Psi_V,\theta_V, \varphi_V)$

are a special set of straight fieldline flux coordinates. To get to straight fieldlines from the `VMEC` system, only the transformation function $\lambda$ is required, and only the poloidal angle has to be transformed. However, one can add an arbitrary periodic function to the transformation and the resulting coordinates are still a straight field line system (see 6.1. "Straight Field-Line Coordinates" in D'haeseleer). The function $G$ is such a in the angle periodic function. As such, it itself can be spectrally decomposed. The ITP internal converter [vmec2boozer](https://gitlab.tugraz.at/plasma/codes/vmec2boozer) or [boozer.py](https://github.com/itpplasma/libneo/blob/main/python/libneo/boozer.py), as well as [`BOOZX`](https://princetonuniversity.github.io/STELLOPT/BOOZ_XFORM) that is part of the [STELLOPT](https://github.com/PrincetonUniversity/STELLOPT) package, all provide this tranformation function as part of their output. There is, however, a bit of a caveat one has to know when using that output.

#### `BOOZX`
The transformation function is given back as spectral coefficients over the new Boozer angles. So the it is not function of the `VMEC` angles ($\theta_V,\varphi_V$) anymore, but of ($\theta_B,\varphi_B$). Comparing to the transformation above, this function is then actually not the transformation 

straight fieldline `VMEC` coordinates &rarr; Boozer coordinates

but rather

straight fieldline `VMEC` coordinates &larr; Boozer coordinates

The idea is to have an easy way to map the toroidal Boozer angle $\varphi_B$ to the geometrical angle $\varphi = \varphi_V$. Therefore `BOOZX` gives the transformation function back as the function $p$

$p \equiv -qG(\Psi_B,\theta_B, \varphi_B)$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Note that $\Psi_B=\Psi_V$).

As a result one can perform the mapping straight-forward

$\varphi = \varphi_B + p(\Psi_B,\theta_B,\varphi_B)$.

#### ITP Boozer

One finds in the `.bc` files of ITP Boozer converters (in contrast to the `.nc` files of `BOOZX`) not the Fourier coefficientss $p_{mn}$ but $v_{mn}$. This $v$ is just $p$ with a normalisation

$v \equiv \frac{n_{\text{per}}}{2\pi}*p(\Psi_B,\theta_B,\varphi_B)$.

The prefactor has to be considered when one again wants to perform the mapping to the geometrical angle as above.

#### The point of it all

The output of the Boozer converter gives back the form of the flux surfaces by providing the modes of the cylindrical coordinates $R$ ($R_{mn}$) and $Z$ ($Z_{mn}$) of the flux surfaces. One can evaluate the spectral series $R(\theta_B,\varphi_B)$ and $Z(\theta_B,\varphi_B)$ then for multiple angles to get points on the surface. To get the surface in real space on needs to go to Cartesian coordinates however, therefore needing the mapping $\varphi_B \rightarrow \varphi$ to calculate

$X = R*\cos{\varphi}$ <br>
$Y = R*\sin{\varphi}$ <br>
$Z = Z$

### Change of spectral modes

#### Inherently different coordinate systems

e.g.
- `VMEC` > `BOOZX`
- `VMEC` > ITP Boozer

One has to recompute the Fourier decomposition as the angles are now different. This requires that one evaluates the spectral series on a grid of points in the angle space i.e. get the values of $f(s,u,v)$. In the direct transformation "`VMEC` > ITP Boozer" done by e.g. [vmec2boozer](https://gitlab.tugraz.at/plasma/codes/vmec2boozer) or [boozer.py](https://github.com/itpplasma/libneo/blob/main/python/libneo/boozer.py) (see python functions of [libneo](https://github.com/itpplasma/libneo/tree/main)), one incorporates the coordinate transfromation into the Fourier integral itself:

$f^B_{mn}(s) = \frac{1}{(2\pi)^2} \int_0^{2\pi}d\theta_V \int_0^{2\pi}d\varphi_V{J_B^*}^{-1} \exp{(-i((m+nq)G)) f(s,\theta_V,\varphi_V)} \exp{(-i(m\theta + n\varphi))}$

This does two things at once:
- Performs the coordinate transformation of the angles ($\theta_V \rightarrow \theta_B, \varphi_V \rightarrow \varphi_B$)
- Changes the Kernel sign (`VMEC` $-$ &rarr; ITP Boozer $+$)

#### Change of Kernel

e.g.
- `BOOZX` > ITP Boozer

When the coordinate system (up to the normalization of the flux label) is the same, and only the sign of the Kernel is the issue, this can simply fixed by making a remapping of the mode coefficients. For example, for converting "`BOOZX` > ITP Boozer" we just

$ f^{\text{ITP Boozer}}_{mn} = f^{\text{BOOZX}}_{m(-n)}$