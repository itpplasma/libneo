# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu January 29 15:17:41 2024

@author: Grassler Georg
"""

# standard modules
import numpy as np
from scipy.interpolate import CubicSpline

class FluxConverter:
    """
    This is a class for the conversion between the poloidal and toroidal flux polflux and torflux of a given flux surface.
    The conversion is based on the safety factor profile q as q = d{toroidal_flux}/d{polodial_flux}.
    Therefore {toroidal_flux} = int_{poloidal_flux_min}^{polodial_flux}q({flux})d{flux}
    """

    def __init__(self, q_profile:np.ndarray, polflux_profile:np.ndarray):
        
        interp_q = CubicSpline(polflux_profile, q_profile, extrapolate=True)
        self.interp_torflux = interp_q.antiderivative()

        # Setup conversion torflux -> polflux using previously made conversion polflux -> torflux
        self.interp_polflux = CubicSpline(self.polflux2torflux(polflux_profile), 
                                        polflux_profile, extrapolate=True) 
        
    def polflux2torflux(self, polflux):
        """
        Converts the poloidal flux polflux to the toroidal flux torflux.
        """

        torflux = self.interp_torflux(polflux)

        return torflux

    def torflux2polflux(self, torflux):
        """
        Converts the toroidal flux torflux to the poloidal flux polflux.
        """

        polflux = self.interp_polflux(torflux)

        return polflux