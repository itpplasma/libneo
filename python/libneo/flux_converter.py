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
    The conversion is based on the safety factor profile q as q = d{Torflux}/d{RibbonPolflux}.
    Therefore {Torflux} = int_{AxisRibbonPolflux}^{RibbonPolflux}q({flux})d{flux} , 
    with {s} = ({RibbonPolflux}-{AxisRibbonPolflux})/({EdgeRibbonPolflux}-{AxisRibbonPolflux})
    gives {Torflux} = ({EdgeRibbonPolflux}-{AxisRibbonPolflux})*int_0^{spol}q(s)ds
          {EdgeTorflux}/({EdgeRibbonPolflux}-{AxisRibbonPolflux}) = int_0^1q(s)ds
    and {stor} = {Tolflux}/{EdgeTorflux} = int_0^{spol}q(s)ds / int_0^1q(s)ds
    
    q_profile: np.ndarray, equidiistant array of the safety factor profile q in terms of poloidal flux

    Furthermore, the absolute fluxes Polflux & Torflux can also be converted with the class using
    the conversion Polflux -> spol -> stor -> Torflux with the methods Polflux2Torflux and Torflux2Polflux.
    For that it requires the absolute poloidal flux as reference,

    AxisPolflux: float, poloidal flux at the magnetic axis (either disk or ribbon flux)
    EdgePolflux: float, poloidal flux at the magnetic edge (either disk or ribbon flux, sames as axis)

    """

    def __init__(self, q_profile:np.ndarray, AxisPolflux:float=np.nan, EdgePolflux:float=np.nan):
        """
        Assumes that the q profile is given on equidistant spol profile from 0 to 1.
        """
        spol_profile = np.linspace(0.0, 1.0, q_profile.shape[0])

        # Setup conversion spol -> stor
        q = CubicSpline(spol_profile, q_profile, extrapolate=True)
        self.__get_Torflux_DeltaRibbonPolflux_ratio = q.antiderivative()
        self.__EdgeTorflux_DeltaRibbonPolflux_ratio = self.__get_Torflux_DeltaRibbonPolflux_ratio(1.0)
        self.__AxisTorflux_DeltaRibbonPolflux_ratio = self.__get_Torflux_DeltaRibbonPolflux_ratio(0.0)

        # Setup conversion stor -> spol using previously made conversion spol -> stor
        self.__spol = CubicSpline(self.spol2stor(spol_profile), 
                                        spol_profile, extrapolate=True)

        # Setup for conversion Polflux <-> Torflux
        self.set_absolute_fluxes(AxisPolflux, EdgePolflux)

    def set_absolute_fluxes(self, AxisPolflux:float, EdgePolflux:float):
        self.__AxisPolflux = AxisPolflux
        self.__EdgePolflux = EdgePolflux
        self.__DeltaPolflux = EdgePolflux - AxisPolflux
        self.__set_EdgeTorflux()

    def __set_EdgeTorflux(self):
        if self.__is_DiskPolflux():
            self.__EdgeTorflux = self.__EdgeTorflux_DeltaRibbonPolflux_ratio * (-self.__DeltaPolflux)
        else:
            self.__EdgeTorflux = self.__EdgeTorflux_DeltaRibbonPolflux_ratio * (+self.__DeltaPolflux)

    def __is_DiskPolflux(self):
        """
        While the RibbonPolflux is zero at the axis, the DiskPolflux is not.
        """
        return self.__AxisPolflux != 0.0

    def spol2stor(self, spol):
        """
        Converts the poloidal flux label spol to the toroidal flux label stor.
        """

        stor = (self.__get_Torflux_DeltaRibbonPolflux_ratio(spol) - self.__AxisTorflux_DeltaRibbonPolflux_ratio)           \
              /(self.__EdgeTorflux_DeltaRibbonPolflux_ratio   - self.__AxisTorflux_DeltaRibbonPolflux_ratio)

        return stor

    def stor2spol(self, stor):
        """
        Converts the toroidal flux label stor to the poloidal flux label spol.
        """

        spol = self.__spol(stor)

        return spol

    def polflux2torflux(self, polflux):
        """
        Converts the poloidal flux polflux to the toroidal flux torflux.
        """

        if np.isnan(self.__AxisPolflux) or np.isnan(self.__EdgePolflux):
            raise ValueError("To use this method, the poloidal flux on axis and edge must be set during \n" + 
                             "initialization or with obj.set_absolute_fluxes(AxisPolflux, EdgePolflux).")

        spol = (polflux-self.__AxisPolflux)/self.__DeltaPolflux
        stor = self.spol2stor(spol)
        torflux = stor * self.__EdgeTorflux

        return torflux

    def torflux2polflux(self, torflux):
        """
        Converts the toroidal flux torflux to the poloidal flux polflux.
        """

        if np.isnan(self.__AxisPolflux) or np.isnan(self.__EdgePolflux):
            raise ValueError("To use this method, the poloidal flux on axis and edge must be set during \n" + 
                             "initialization or with obj.set_absolute_fluxes(AxisPolflux, EdgePolflux).")

        stor = torflux / self.__EdgeTorflux
        spol = self.stor2spol(stor)
        polflux = spol * self.__DeltaPolflux + self.__AxisPolflux

        return polflux