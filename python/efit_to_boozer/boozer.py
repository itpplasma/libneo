#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:17:41 2016

@author: Christopher Albert
"""
import numpy as np
import _efit_to_boozer

length_cgs_to_si = 1e-2
debug = False


def init_efit_to_boozer(construction=0, **kwargs):
    """Initialize native magdata, raising on invalid input or equilibrium."""
    try:
        _efit_to_boozer.efit_to_boozer.init(construction=construction, **kwargs)
    except ValueError as error:
        raise RuntimeError(f"efit_to_boozer init failed: {error}") from error


def _construction_options(construction, method, gfile):
    """Preserve the numeric selector and optional native deck configuration."""
    from pathlib import Path
    import shlex

    if construction not in (0, 1, 2):
        raise ValueError("EFIT construction must be 0 (deck), 1 (contour), or 2 (legacy)")
    selected, options = construction or 1, {}
    deck = Path("efit_to_boozer.inp")
    if construction == 0 and gfile is None and deck.exists():
        records = deck.read_text().splitlines()
        if len(records) > 7:
            selected = int(records[7].split()[0])
        if len(records) > 8:
            options["route"] = int(records[8].split()[0])
        if len(records) > 9:
            options["allow_poor_equilibrium"] = records[9].split()[0].lower() in ("t", ".true.", "true")
        if len(records) > 10:
            options["override_reason"] = shlex.split(records[10], comments=True)[0]
        if len(records) > 11:
            options["n_theta"] = int(records[11].split()[0])
        if len(records) > 13:
            options["exp_bpol"] = int(records[13].split()[0])
    if selected not in (1, 2):
        raise ValueError("EFIT deck construction must be 1 (contour) or 2 (legacy)")
    if method is None:
        method = "legacy" if selected == 2 else "contour"
    if method not in ("contour", "legacy"):
        raise ValueError("EFIT method must be contour or legacy")
    return method, options


def _contour_input(gfile, psimax):
    from pathlib import Path
    import shlex

    from_deck = gfile is None
    if from_deck:
        records = Path("field_divB0.inp").read_text().splitlines()
        if len(records) < 7:
            raise ValueError("field_divB0.inp must contain the g-file on record 7")
        gfile = shlex.split(records[6], comments=True)[0]
    if from_deck and psimax is None and Path("efit_to_boozer.inp").exists():
        records = Path("efit_to_boozer.inp").read_text().splitlines()
        if len(records) < 7:
            raise ValueError("efit_to_boozer.inp must contain seven records")
        boundary = float(records[6].split()[0].replace("D", "E").replace("d", "e"))
        if not np.isfinite(boundary):
            raise ValueError("psimax must be finite")
        if boundary != 1.0e10 and boundary < 1.0e299:
            psimax = boundary
    return str(Path(gfile).resolve()), psimax


def get_contour_transform(gfile, s, mpol, *, n_theta=1024, n_q=128, route=1,
                          psimax=None, allow_poor_equilibrium=False,
                          override_reason="", exp_bpol=-1):
    """Return native contour spectra and flux functions in counterclockwise theta.

    Lengths in ``spectra`` are cm, fields G and currents G cm. ``axis`` is
    in metres. The v coefficients obey phi_B = phi - 2*pi*v. ``psitor_max``
    retains the native oriented poloidal-flux sign, in G cm^2 per radian.
    The five route numbers are those of the Fortran ``efit_sfl`` module.
    ``exp_bpol`` is -1 to infer units, 0 for Wb/rad, or 1 for Wb; an
    unidentified source requires explicit units even with a quality override.
    """
    return _get_contour_transform(gfile, s, mpol, n_theta=n_theta, n_q=n_q,
                                  route=route, psimax=psimax,
                                  allow_poor_equilibrium=allow_poor_equilibrium,
                                  override_reason=override_reason, exp_bpol=exp_bpol)


def get_contour_transform_extended(gfile, s, mpol, **options):
    """Add native pressure derivative, volume factor, radius and orientation.

    ``pp_file`` and ``sqrtg00`` retain native Gaussian units and signs;
    ``rminor`` is in cm and ``sig_theta`` is the native poloidal orientation.
    Other entries and options have the contract of :func:`get_contour_transform`.
    """
    return _get_contour_transform(gfile, s, mpol, extended_metadata=True, **options)


def _get_contour_transform(gfile, s, mpol, *, n_theta=1024, n_q=128, route=1,
                          psimax=None, allow_poor_equilibrium=False,
                          override_reason="", exp_bpol=-1, extended_metadata=False):
    """Read either flat native layout and apply the shared CCW orientation."""
    gfile, psimax = _contour_input(gfile, psimax)
    s = np.asarray(s, dtype=float)
    if s.ndim != 1 or not len(s) or not np.all(np.isfinite(s)):
        raise ValueError("s must be a nonempty finite one-dimensional array")
    if psimax is not None and not np.isfinite(psimax):
        raise ValueError("psimax must be finite")
    if allow_poor_equilibrium and not override_reason.strip():
        raise ValueError("allow_poor_equilibrium requires an explicit override_reason")
    native = _efit_to_boozer.efit_to_boozer
    transform = (native.contour_transform_extended if extended_metadata
                 else native.contour_transform)
    result = transform(gfile, s, mpol, n_theta, n_q, route, psimax is not None,
                       0.0 if psimax is None else psimax, allow_poor_equilibrium,
                       override_reason, exp_bpol=exp_bpol)
    if extended_metadata:
        coefficients, fluxes, state, physical, rminor, code, overridden, message = result
    else:
        coefficients, fluxes, state, code, overridden, message = result
    if isinstance(message, bytes):
        message = message.decode().strip()
    else:
        message = str(message).strip()
    if code and not overridden:
        raise ValueError(message or f"EFIT contour transform failed with status {code}")
    # Native theta follows poloidal motion at increasing phi. Reorient its
    # sine modes and poloidal covariant component together for this Python API.
    coefficients = np.array(coefficients, copy=True)
    coefficients[[1, 3, 5, 7], :, :] *= state[3]
    metadata = dict(
        axis=np.asarray(state[:2])*length_cgs_to_si,
        iota=np.asarray(fluxes[0])*state[3],
        curr_i=np.asarray(fluxes[1])*state[3],
        curr_g=np.asarray(fluxes[2]),
        psitor_max=float(state[2]),
        spectra=coefficients,
        status_code=int(code), overridden=bool(overridden), message=message,
        gfile=gfile,
    )
    if extended_metadata:
        metadata.update(pp_file=np.asarray(physical[0]), sqrtg00=np.asarray(physical[1]),
                        rminor=float(rminor), sig_theta=float(state[3]))
    return metadata


def get_boozer_transform(stor, num_theta, construction=0, *, method=None, gfile=None,
                         n_theta=None, n_q=128, psimax=None,
                         allow_poor_equilibrium=None, override_reason=None,
                         return_metadata=False, route=None, exp_bpol=None):
    """Return periodic Boozer geometry splines, using contour construction.

    The seven-item return tuple and the unused final ``stor`` sentinel are
    retained. ``method='legacy'`` selects historical field-line integration.
    ``return_metadata=True`` additionally returns the contour flux functions.
    The contour grid has at least four points per retained Fourier mode.
    """
    from scipy.interpolate import CubicSpline

    method, options = _construction_options(construction, method, gfile)
    for key, value in dict(n_theta=n_theta, route=route, exp_bpol=exp_bpol,
                           allow_poor_equilibrium=allow_poor_equilibrium,
                           override_reason=override_reason).items():
        if value is not None:
            options[key] = value
    options.setdefault("n_theta", 1024)
    options["n_theta"] = max(options["n_theta"], 4*num_theta)
    if method == "legacy":
        if return_metadata:
            raise ValueError("Contour metadata is unavailable for method='legacy'")
        return _get_legacy_boozer_transform(np.asarray(stor), num_theta)
    if method != "contour":
        raise ValueError("EFIT method must be contour or legacy")
    stor = np.asarray(stor, dtype=float)
    if stor.ndim != 1 or len(stor) < 2:
        raise ValueError("stor must contain surfaces and a final sentinel")
    if not isinstance(num_theta, (int, np.integer)) or num_theta < 4:
        raise ValueError("num_theta must be an integer >= 4")
    data = get_contour_transform(
        gfile, stor[:-1], num_theta, n_q=n_q, psimax=psimax, **options,
    )
    theta = np.linspace(0.0, 2*np.pi, 2*num_theta + 1)
    modes = np.arange(num_theta + 1)[:, None]
    co, si = np.cos(modes*theta), np.sin(modes*theta)
    coefficients = data["spectra"]
    samples = np.einsum("fms,mt->sft", coefficients[0::2], co)
    samples += np.einsum("fms,mt->sft", coefficients[1::2], si)
    r_of_thb, dth_of_thb, G_of_thb, B0 = [], [], [], []
    theta_geoms = np.zeros((len(stor)-1, len(theta)))
    theta_symflux = theta[:-1]
    for k, (R, Z, v, B) in enumerate(samples):
        R = R*length_cgs_to_si - data["axis"][0]
        Z = Z*length_cgs_to_si - data["axis"][1]
        geom = np.unwrap(np.arctan2(Z, R))
        radius = np.hypot(R, Z)
        dtheta = geom - theta
        nu = -2*np.pi*v
        for values in (radius, dtheta, nu, B):
            values[-1] = values[0]
        r_of_thb.append(CubicSpline(theta, radius, bc_type="periodic"))
        dth_of_thb.append(CubicSpline(theta, dtheta, bc_type="periodic"))
        G_of_thb.append(CubicSpline(theta, nu, bc_type="periodic"))
        B0.append(CubicSpline(theta, B, bc_type="periodic"))
        theta_geoms[k] = theta + dtheta
        theta_symflux = theta[:-1] - data["iota"][k]*nu[:-1]
    result = (r_of_thb, dth_of_thb, G_of_thb,
              np.broadcast_to(theta, theta_geoms.shape).copy(), theta_geoms,
              theta_symflux, B0)
    return (result, data) if return_metadata else result


def get_angles_G_and_B0(s, num_theta, construction=0, *, method=None, **kwargs):
    method, _ = _construction_options(construction, method, kwargs.get("gfile"))
    if method == "legacy":
        return _get_legacy_angles_G_and_B0(s, num_theta)
    result = get_boozer_transform(np.array([s, 2.0]), num_theta,
                                  construction, method=method, **kwargs)
    radius, dtheta, nu, theta, geom, symflux, B = result
    return (radius[0](theta[0]), theta[0], geom[0], symflux,
            nu[0](theta[0]), B[0](theta[0]))


def get_B0_of_s_theta_boozer(stor, num_theta, construction=0, *, method=None, **kwargs):
    return get_boozer_transform(stor, num_theta, construction, method=method, **kwargs)[-1]


def get_magnetic_axis(construction=0, *, method=None, gfile=None, **kwargs):
    method, options = _construction_options(construction, method, gfile)
    options.update(kwargs)
    if method == "legacy":
        return _get_legacy_magnetic_axis()
    if method != "contour":
        raise ValueError("EFIT method must be contour or legacy")
    data = get_contour_transform(gfile, [0.5], 0, **options)
    return tuple(data["axis"])


def _get_legacy_boozer_transform(stor, num_theta):
    from scipy.interpolate import CubicSpline

    _efit_to_boozer.efit_to_boozer.init(method='legacy')

    ns = stor.shape[0]

    G_of_thb = []
    dth_of_thb = []
    r_of_thb = []
    B0 = []
    theta_boozers = np.zeros((ns - 1, 2*num_theta+1))
    theta_geoms = np.zeros((ns - 1, 2*num_theta+1))
    for ks in np.arange(ns - 1):
        r, th_boozers, th_geoms, theta_symflux, G_s, B0_temp = _get_legacy_angles_G_and_B0(stor[ks], num_theta)
        dth = th_geoms - th_boozers

        G_of_thb.append(CubicSpline(th_boozers, G_s, bc_type="periodic"))
        dth_of_thb.append(CubicSpline(th_boozers, dth, bc_type="periodic"))
        r_of_thb.append(CubicSpline(th_boozers, r, bc_type="periodic"))
        B0.append(CubicSpline(th_boozers, B0_temp, bc_type='periodic'))
        theta_boozers[ks, :] = th_boozers
        theta_geoms[ks,:] = th_geoms

    return r_of_thb, dth_of_thb, G_of_thb, theta_boozers, theta_geoms, theta_symflux, B0


def _get_legacy_angles_G_and_B0(s, num_theta):
    psi = np.array(0.0)
    R_axis, Z_axis = get_magnetic_axis(method='legacy')

    th_boozers = []
    th_geoms = []
    r = []
    Gs = []
    B0mod = []
    theta_symflux = 2 * np.pi * np.arange(2 * num_theta) / (2 * num_theta)

    for th_symflux in theta_symflux:
        inp_label = 1
        (q, _, _, _, _,
        B0mod_temp, _,
        R, _, _,
        Z, _, _,
        G, _, _,
        ) = _efit_to_boozer.efit_to_boozer.magdata(inp_label, s, psi, th_symflux)
        Gs.append(G)
        th_boozers.append(th_symflux + G / q)
        th_geoms.append(np.arctan2(
          Z*length_cgs_to_si - Z_axis,
          R*length_cgs_to_si - R_axis))
        r.append(np.sqrt((R*length_cgs_to_si - R_axis)**2 + (Z*length_cgs_to_si - Z_axis)**2))
        B0mod.append(B0mod_temp)

    Gs = np.array(Gs)

    B0mod = np.array(B0mod)

    r = np.array(r)

    th_boozers = np.unwrap(th_boozers)
    monotone_sign_boozers = np.sign(np.mean(np.sign(np.diff(th_boozers))))

    th_geoms = np.unwrap(th_geoms)
    monotone_sign_geoms = np.sign(np.mean(np.sign(np.diff(th_geoms))))

    if monotone_sign_boozers != monotone_sign_geoms:
        th_boozers = np.unwrap(-th_boozers)

    if monotone_sign_geoms < 0:
        r = r[::-1]
        th_boozers = np.unwrap(th_boozers[::-1])
        th_geoms = np.unwrap(th_geoms[::-1])
        Gs = Gs[::-1]
        B0mod = B0mod[::-1]

    # Reverse samples and their abscissae together before closing the period.
    # Closing only the values first shifts reversed data by one poloidal sample.
    r = np.append(r, r[0])
    Gs = np.append(Gs, Gs[0])
    B0mod = np.append(B0mod, B0mod[0])
    th_boozers = np.append(th_boozers, th_boozers[0] + 2 * np.pi)
    th_geoms = np.append(th_geoms, th_geoms[0] + 2 * np.pi)

    return r, th_boozers, th_geoms, theta_symflux, Gs, B0mod


def _get_legacy_B0_of_s_theta_boozer(stor, num_theta):
  from scipy.interpolate import CubicSpline

  _efit_to_boozer.efit_to_boozer.init(method='legacy')

  ns = stor.shape[0]

  B0 = []
  for ks in np.arange(ns - 1):
    psi = np.array(0.0)
    s = np.array(stor[ks])

    th_boozers = []
    B0mod = []
    for th_symflux in 2 * np.pi * np.arange(2 * num_theta) / (2 * num_theta):
      inp_label = 1
      (
        q,
        _,
        _,
        _,
        _,
        B0mod_temp,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        G,
        _,
        _,
      ) = _efit_to_boozer.efit_to_boozer.magdata(
        inp_label, s, psi, th_symflux)
      th_boozers.append(th_symflux + G / q)
      B0mod.append(B0mod_temp)

    th_boozers = np.unwrap(th_boozers)
    if not all(np.diff(th_boozers)>0):
      raise ValueError("Boozer angles are not monotonically increasing")
    th_boozers = np.append(th_boozers, th_boozers[0] + 2 * np.pi)

    B0mod = np.array(B0mod)
    B0mod = np.append(B0mod, B0mod[0])
    B0.append(CubicSpline(th_boozers, B0mod, bc_type='periodic'))

  return B0


def get_boozer_harmonics(f, stor, num_theta, num_phi, m0b, n, dth_of_thb, G_of_thb):
  """
  f(stor, th_geom, ph) of normalized toroidal flux and geometric angles

  Returns: Fourier harmonics fmn in terms of Boozer angles
  Set num_phi=1 for axisymmetric case.
  """

  ns = stor.shape[0]

  fmn = np.zeros((ns - 1, 2 * m0b + 1), dtype=complex)

  th_geoms, phs = get_thgeoms_phs(ns, num_theta, num_phi, dth_of_thb, G_of_thb)


  for kth in np.arange(num_theta):
    if debug: print(f"kth = {kth}/{num_theta-1}")
    thb = kth * 2 * np.pi / (num_theta)
    for kph in np.arange(num_phi):
      phb = kph * 2 * np.pi / (num_phi)
      f_values = f(stor[:-1], th_geoms[:, kth], phs[:, kth, kph])

      for m in np.arange(-m0b, m0b + 1):
        # Take a sum over all theta values here
        fmn[:, m + m0b] += f_values * np.exp(-1j * (m * thb + n * phb))

  return fmn / ((num_theta) * (num_phi))

def get_thgeoms_phs(ns, num_theta, num_phi, dth_of_thb, G_of_thb):

  th_geoms = np.zeros((ns - 1, num_theta))
  phs = np.zeros((ns - 1, num_theta, num_phi))

  # Loop over flux surface index
  for ks in np.arange(ns - 1):
    if debug: print(f"ks = {ks}/{ns-2}")
    for kth in np.arange(num_theta):
      thb = kth * 2 * np.pi / num_theta
      th_geoms[ks, kth] = thb + dth_of_thb[ks](thb)
      G = G_of_thb[ks](thb)
      for kph in np.arange(num_phi):
        phb = kph * 2 * np.pi / num_phi
        phs[ks, kth, kph] = phb - G

  return th_geoms, phs

def get_boozer_harmonics_divide_f_by_B0(f, stor, num_theta, num_phi, m0b, n, dth_of_thb, G_of_thb, **kwargs):
  """
  f(stor, th_geom, ph) of normalized toroidal flux and geometric angles

  Returns: Fourier harmonics fmn in terms of Boozer angles, i.e. fmn(s,m) for a given
  toroidal mode number n
  """

  ns = stor.shape[0]
  fmn = np.zeros((ns - 1, 2 * m0b + 1), dtype=complex)

  th_geoms, phs = get_thgeoms_phs(ns, num_theta, num_phi, dth_of_thb, G_of_thb)
  B0 = get_B0_of_s_theta_boozer(stor, num_theta, **kwargs)

  for kth in np.arange(num_theta):
    if debug: print(f"kth = {kth}/{num_theta-1}")
    thb = kth * 2 * np.pi / num_theta
    B0_arr = np.array([spline(thb) for spline in B0])

    for kph in np.arange(num_phi):
      phb = kph * 2 * np.pi / num_phi

      f_values = np.array(f(stor[:-1], th_geoms[:, kth], phs[:, kth, kph]), dtype=complex) \
        / B0_arr

      for m in np.arange(-m0b, m0b + 1):
        # Take a sum over all theta values here
        fmn[:, m + m0b] += f_values * np.exp(-1j * (m * thb + n * phb))

  return fmn / (num_theta * num_phi)

def get_boozer_harmonics_divide_f_by_B0_1D(f, stor, num_theta, num_phi, m0b, n, dth_of_thb, G_of_thb, **kwargs):
  """
  f_n(stor, th_geom) of normalized toroidal flux and geometric poloidal angle for given
  toroidal mode number n (also in geometric angle). Takes also the transformation from
  geometric to Boozer angle in toroidal direction into account (additional phase factor
  depending on G(\vatheta_B)).

  Returns: Fourier harmonics fmn in terms of Boozer angles, i.e. fmn(s,m) for a given
  toroidal mode number n
  """

  ns = stor.shape[0]
  fmn = np.zeros((ns - 1, 2 * m0b + 1), dtype=complex)

  th_geoms, phs = get_thgeoms_phs(ns, num_theta, num_phi, dth_of_thb, G_of_thb)
  B0 = get_B0_of_s_theta_boozer(stor, num_theta, **kwargs)

  for kth in np.arange(num_theta):
    if debug: print(f"kth = {kth}/{num_theta-1}")
    thb = kth * 2 * np.pi / num_theta
    B0_arr = np.array([spline(thb) for spline in B0])

    G = np.array([G_of_thb[ks](thb) for ks in np.arange(ns - 1)])

    f_values = np.array(f(stor[:-1], th_geoms[:, kth], phs[:, kth, 0]), dtype=complex) \
      / B0_arr * np.exp(-1j * n * G)

    for m in np.arange(-m0b, m0b + 1):
      # Take a sum over all theta values here
      fmn[:, m + m0b] += f_values * np.exp(-1j * (m * thb))

  return fmn / (num_theta * num_phi)

def _get_legacy_magnetic_axis():
  psi = np.array(0.0)
  theta = np.array(0.0)
  si = np.array(0.0)

  inp_label = 1
  (
    _,
    _,
    _,
    _,
    _,
    _,
    _,
    R,
    _,
    _,
    Z,
    _,
    _,
    _,
    _,
    _,
  ) = _efit_to_boozer.efit_to_boozer.magdata(inp_label, si, psi, theta)

  R_axis = R * length_cgs_to_si
  Z_axis = Z * length_cgs_to_si

  return R_axis, Z_axis
