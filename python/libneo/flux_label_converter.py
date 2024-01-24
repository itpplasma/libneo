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
    Therefore {toroidal_flux} = int_-inf^{polodial_flux}q({flux})d{flux} - , with {s} = {flux}/{flux_max}
    gives {toroidal_flux} = {flux_max}*int_-inf^{s_pol}q(s)ds
    and {s_tor} = {todroidal_flux}/{toroidal_flux_max} = int_-inf^{s_pol}q(s)ds / int_-inf^1q(s)ds
    """

    def __init__(self, q_profile:np.ndarray):
        """
        Assumes that the q profile is given on equidistant s_pol profile from 0 to 1.
        """
        spol_list = np.linspace(0.0, 1.0, q_profile.shape[0])

        # Setup conversion spol -> stor
        interp_q = CubicSpline(spol_list, q_profile, extrapolate=True)
        self.interp_psitor = interp_q.antiderivative()
        self.psitor_max = self.interp_psitor(1.0)
        self.psitor_min = self.interp_psitor(0.0)

        # Setup conversion stor -> spol using previously made conversion spol -> stor
        self.interp_spol = CubicSpline(self.spol2stor(spol_list), 
                                        spol_list, extrapolate=True)

    def spol2stor(self, spol):
        """
        Converts the poloidal flux label s_pol to the toroidal flux label s_tor.
        """

        stor = (self.interp_psitor(spol)-self.psitor_min)/(self.psitor_max - self.psitor_min)

        return stor

    def stor2spol(self, stor):
        """
        Converts the toroidal flux label s_tor to the poloidal flux label s_pol.
        """

        spol = self.interp_spol(stor)

        return spol

################################################################################
################################################################################
################################################################################

if __name__=='__main__':
    print('This file contains a class for the conversion between the poloidal and') 
    print('toroidal flux label s_pol and s_tor of a given flux surface.')