"""
Created on Thu January 18 15:17:41 2024

@author: Grassler Georg
"""

import numpy as np
from scipy.interpolate import CubicSpline


class FluxConverter:
    """
    This is a class for the conversion between the poloidal and toroidal flux of a given flux surface.
    """

    def __init__(self, q_profile: np.ndarray, spol_profile: np.ndarray = np.nan,
                                              axis_polflux: float = np.nan, 
                                              edge_polflux: float = np.nan):
        """
        q_profile: np.ndarray, safety factor profile q in terms of poloidal label spol
        spol_profile: (optional) np.ndarry, if not provided, an equidistant profile from 0 to 1 is assumed
        axis_polflux: (optional) float, poloidal flux at the magnetic axis (either disk or ribbon flux)
        edge_polflux: (optional) float, poloidal flux at the magnetic edge (either disk or ribbon flux, sames as axis)
        """
        if np.any(np.isnan(spol_profile)):
            spol_profile = np.linspace(0.0, 1.0, q_profile.shape[0])

        # Setup conversion spol -> stor
        q = CubicSpline(spol_profile, q_profile, extrapolate=True)
        self.__get_torflux_delta_ribbon_polflux_ratio = q.antiderivative()
        self.__edge_torflux_delta_ribbon_polflux_ratio = self.__get_torflux_delta_ribbon_polflux_ratio(1.0)
        self.__axis_torflux_delta_ribbon_polflux_ratio = self.__get_torflux_delta_ribbon_polflux_ratio(0.0)

        # Setup conversion stor -> spol using previously made conversion spol -> stor
        self.__spol = CubicSpline(self.spol2stor(spol_profile),
                                  spol_profile, extrapolate=True)

        # Setup for conversion Polflux <-> Torflux
        self.set_absolute_fluxes(axis_polflux, edge_polflux)

    def set_absolute_fluxes(self, axis_polflux: float, edge_polflux: float):
        self.__axis_polflux = axis_polflux
        self.__edge_polflux = edge_polflux
        self.__delta_polflux = edge_polflux - axis_polflux
        self.__set_edge_torflux()

    def __set_edge_torflux(self):
        if self.__is_disk_polflux():
            self.__edge_torflux = self.__edge_torflux_delta_ribbon_polflux_ratio * (-self.__delta_polflux)
        else:
            self.__edge_torflux = self.__edge_torflux_delta_ribbon_polflux_ratio * (+self.__delta_polflux)

    def __is_disk_polflux(self):
        """
        While the RibbonPolflux is zero at the axis, the disk_polflux is not.
        """
        return self.__axis_polflux != 0.0

    def spol2stor(self, spol):

        stor = (self.__get_torflux_delta_ribbon_polflux_ratio(spol) - self.__axis_torflux_delta_ribbon_polflux_ratio) \
               / (self.__edge_torflux_delta_ribbon_polflux_ratio - self.__axis_torflux_delta_ribbon_polflux_ratio)

        return stor

    def stor2spol(self, stor):

        spol = self.__spol(stor)

        return spol

    def stor2sqrtspol(self, stor):

        spol = self.stor2spol(stor)
        sqrtspol = np.sqrt(spol)

        return sqrtspol

    def sqrtspol2stor(self, sqrtspol):

        spol = sqrtspol**2
        stor = self.spol2stor(spol)

        return stor

    def polflux2torflux(self, polflux):

        if np.isnan(self.__axis_polflux) or np.isnan(self.__edge_polflux):
            raise ValueError(
                "To use this method, the poloidal flux on axis and edge must be set during \n" +
                "initialization or with obj.set_absolute_fluxes(axis_polflux, edge_polflux).")

        spol = (polflux - self.__axis_polflux) / self.__delta_polflux
        stor = self.spol2stor(spol)
        torflux = stor * self.__edge_torflux

        return torflux

    def torflux2polflux(self, torflux):

        if np.isnan(self.__axis_polflux) or np.isnan(self.__edge_polflux):
            raise ValueError(
                "To use this method, the poloidal flux on axis and edge must be set during \n" +
                "initialization or with obj.set_absolute_fluxes(axis_polflux, edge_polflux).")

        stor = torflux / self.__edge_torflux
        spol = self.stor2spol(stor)
        polflux = spol * self.__delta_polflux + self.__axis_polflux

        return polflux
