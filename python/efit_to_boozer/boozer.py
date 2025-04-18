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


def get_boozer_transform(stor, num_theta):
    from scipy.interpolate import CubicSpline

    _efit_to_boozer.efit_to_boozer.init()

    ns = stor.shape[0]

    G_of_thb = []
    dth_of_thb = []
    r_of_thb = []
    B0 = []
    theta_boozers = np.zeros((ns - 1, 2*num_theta+1))
    theta_geoms = np.zeros((ns - 1, 2*num_theta+1))
    for ks in np.arange(ns - 1):
        r, th_boozers, th_geoms, theta_symflux, G_s, B0_temp = get_angles_G_and_B0(stor[ks], num_theta)
        sign, dth, th_boozers, th_geoms = get_sign_dependent_thetas(th_geoms, th_boozers)

        G_of_thb.append(CubicSpline(th_boozers, G_s, bc_type="periodic"))
        dth_of_thb.append(CubicSpline(th_boozers, dth, bc_type="periodic"))
        r_of_thb.append(CubicSpline(th_boozers, r, bc_type="periodic"))
        B0.append(CubicSpline(th_boozers, B0_temp, bc_type='periodic'))
        theta_boozers[ks, :] = th_boozers
        theta_geoms[ks,:] = th_geoms

    return sign, r_of_thb, dth_of_thb, G_of_thb, theta_boozers, theta_geoms, theta_symflux, B0


def get_angles_G_and_B0(s, num_theta):
    psi = np.array(0.0)
    R_axis, Z_axis = get_magnetic_axis()

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
    Gs = np.append(Gs, Gs[0])

    B0mod = np.array(B0mod)
    B0mod = np.append(B0mod, B0mod[0])

    r = np.array(r)
    r = np.append(r, r[0])

    return r, th_boozers, th_geoms, theta_symflux, Gs, B0mod


def get_sign_dependent_thetas(th_geoms, th_boozers):

  th_boozers = np.unwrap(th_boozers)
  monotone_sign_boozers = np.sign(np.mean(np.sign(np.diff(th_boozers))))

  th_geoms = np.unwrap(th_geoms)
  monotone_sign_geoms = np.sign(np.mean(np.sign(np.diff(th_geoms))))

  th_boozers = np.append(th_boozers, th_boozers[0] + monotone_sign_boozers * 2 * np.pi)
  th_geoms = np.append(th_geoms, th_geoms[0] + monotone_sign_geoms * 2 * np.pi)

  if monotone_sign_boozers == monotone_sign_geoms:
    sign = -1
  else:
    sign = 1

  dth = th_geoms + sign * th_boozers

  return sign, dth, th_boozers, th_geoms



def get_B0_of_s_theta_boozer(stor, num_theta):
  from scipy.interpolate import CubicSpline

  _efit_to_boozer.efit_to_boozer.init()

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
    monotone_sign_boozers = np.sign(np.mean(np.sign(np.diff(th_boozers))))
    th_boozers = np.append(th_boozers, th_boozers[0] + monotone_sign_boozers * 2 * np.pi)

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

def get_boozer_harmonics_divide_f_by_B0(f, stor, num_theta, num_phi, m0b, n, dth_of_thb, G_of_thb):
  """
  f(stor, th_geom, ph) of normalized toroidal flux and geometric angles

  Returns: Fourier harmonics fmn in terms of Boozer angles, i.e. fmn(s,m) for a given
  toroidal mode number n
  """

  ns = stor.shape[0]
  fmn = np.zeros((ns - 1, 2 * m0b + 1), dtype=complex)

  th_geoms, phs = get_thgeoms_phs(ns, num_theta, num_phi, dth_of_thb, G_of_thb)
  B0 = get_B0_of_s_theta_boozer(stor, num_theta)

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

def get_boozer_harmonics_divide_f_by_B0_1D(f, stor, num_theta, num_phi, m0b, n, dth_of_thb, G_of_thb):
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
  B0 = get_B0_of_s_theta_boozer(stor, num_theta)

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

def get_magnetic_axis():
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
