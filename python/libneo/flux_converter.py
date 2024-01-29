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
    This is a class for the conversion between the poloidal and toroidal flux label psipol and phitor of a given flux surface.
    The conversion is based on the safety factor profile q as q = d{toroidal_flux}/d{polodial_flux}.
    Therefore {toroidal_flux} = int_{poloidal_flux_min}^{polodial_flux}q({flux})d{flux}
    """

    def __init__(self, q_profile:np.ndarray, psipol_profile:np.ndarray):
        
        interp_q = CubicSpline(psipol_profile, q_profile, extrapolate=True)
        self.interp_phitor = interp_q.antiderivative()

        # Setup conversion phitor -> psipol using previously made conversion psipol -> phitor
        self.interp_psipol = CubicSpline(self.psipol2phitor(psipol_profile), 
                                        psipol_profile, extrapolate=True) 
        
    def psipol2phitor(self, psipol):
        """
        Converts the poloidal flux psipol to the toroidal flux phitor.
        """

        phitor = self.interp_phitor(psipol)

        return phitor

    def phitor2psipol(self, phitor):
        """
        Converts the toroidal flux phitor to the poloidal flux psipol.
        """

        psipol = self.interp_psipol(phitor)

        return psipol

################################################################################
################################################################################
################################################################################

if __name__=='__main__':
    print('This file contains a class for the conversion between the poloidal and') 
    print('toroidal flux psipol and phitor of a given flux surface.')