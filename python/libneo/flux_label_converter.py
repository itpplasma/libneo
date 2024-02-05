# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu January 18 15:17:41 2024

@author: Grassler Georg
"""

import numpy as np
from scipy.interpolate import CubicSpline

class FluxLabelConverter:
    """
    This is a class for the conversion between the poloidal and toroidal flux label s_pol and s_tor of a given flux surface.
    The conversion is based on the safety factor profile q as q = d{toroidal_flux}/d{polodial_flux}.
    Therefore {toroidal_flux} = int_{poloidal_flux_min}^{polodial_flux}q({flux})d{flux} , with {s} = ({poloidal_flux}-{poloidal_flux_min})/({poloidal_flux_max}-{poloidal_flux_min})
    gives {toroidal_flux} = ({poloidal_flux_max}-{poloidal_flux_min})*int_0^{s_pol}q(s)ds
          {toroidal_flux_max}/({poloidal_flux_max}-{poloidal_flux_min}) = int_0^1q(s)ds
    and {s_tor} = {todroidal_flux}/{toroidal_flux_max} = int_0^{s_pol}q(s)ds / int_0^1q(s)ds
    """

    def __init__(self, q_profile:np.ndarray):
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