# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu January 18 15:17:41 2024

@author: Grassler Georg
"""

import numpy as np
from scipy.interpolate import CubicSpline

class FluxConverter:
    """
    This is a class for the conversion between the poloidal and toroidal flux of a given flux surface.
    The conversion is based on the safety factor profile q as q = d{toroidal_flux}/d{polodial_flux}.
    Therefore {toroidal_flux} = int_{poloidal_flux_min}^{polodial_flux}q({flux})d{flux} , with {s} = ({poloidal_flux}-{poloidal_flux_min})/({poloidal_flux_max}-{poloidal_flux_min})
    gives {toroidal_flux} = ({poloidal_flux_max}-{poloidal_flux_min})*int_0^{s_pol}q(s)ds
          {toroidal_flux_max}/({poloidal_flux_max}-{poloidal_flux_min}) = int_0^1q(s)ds
    and {s_tor} = {todroidal_flux}/{toroidal_flux_max} = int_0^{s_pol}q(s)ds / int_0^1q(s)ds
    
    q_profile: np.ndarray, equidiistant array of the safety factor profile q in terms of poloidal flux

    Furthermore, the absolute fluxes polflux & torfluc can also be converted with the class using
    the conversion polflux -> s_pol -> s_tor -> torflux with the methods polflux2torflux and torflux2polflux.
    For that it requires the absolute poloidal flux as reference,

    axis_polflux: float, poloidal flux at the magnetic axis
    edge_polflux: float, poloidal flux at the magnetic edge

    """

    def __init__(self, q_profile:np.ndarray, axis_polflux:float=np.nan, edge_polflux:float=np.nan):
        """
        Assumes that the q profile is given on equidistant s_pol profile from 0 to 1.
        """
        spol_list = np.linspace(0.0, 1.0, q_profile.shape[0])

        # Setup conversion spol -> stor
        interp_q = CubicSpline(spol_list, q_profile, extrapolate=True)
        self.interp_torflux_over_delta_polflux = interp_q.antiderivative()
        self.torflux_max_over_delta_polflux = self.interp_torflux_over_delta_polflux(1.0)
        self.torflux_min_over_delta_polflux = self.interp_torflux_over_delta_polflux(0.0)

        # Setup conversion stor -> spol using previously made conversion spol -> stor
        self.interp_spol = CubicSpline(self.spol2stor(spol_list), 
                                        spol_list, extrapolate=True)

        # Setup for conversion polflux <-> torflux
        self.axis_polflux = axis_polflux
        self.edge_polflux = edge_polflux
        self.delta_polflux = edge_polflux - axis_polflux
        self.torflux_max = self.torflux_max_over_delta_polflux * self.delta_polflux

    def spol2stor(self, spol):
        """
        Converts the poloidal flux label s_pol to the toroidal flux label s_tor.
        """

        stor = (self.interp_torflux_over_delta_polflux(spol)-self.torflux_min_over_delta_polflux)/(self.torflux_max_over_delta_polflux - self.torflux_min_over_delta_polflux)

        return stor

    def stor2spol(self, stor):
        """
        Converts the toroidal flux label s_tor to the poloidal flux label s_pol.
        """

        spol = self.interp_spol(stor)

        return spol

    def polflux2torflux(self, polflux):
        """
        Converts the poloidal flux polflux to the toroidal flux torflux.
        """

        if np.isnan(self.axis_polflux) or np.isnan(self.edge_polflux):
            raise ValueError("The poloidal flux axis and edge must be set in initialization to use this method.")

        spol = (polflux-self.axis_polflux)/self.delta_polflux
        stor = self.spol2stor(spol)
        torflux = stor * self.torflux_max

        return torflux

    def torflux2polflux(self, torflux):
        """
        Converts the toroidal flux torflux to the poloidal flux polflux.
        """

        if np.isnan(self.axis_polflux) or np.isnan(self.edge_polflux):
            raise ValueError("The poloidal flux axis and edge must be set in initialization to use this method.")

        stor = torflux / self.torflux_max
        spol = self.stor2spol(stor)
        polflux = spol * self.delta_polflux + self.axis_polflux

        return polflux