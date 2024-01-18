# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu January 18 15:17:41 2024

@author: Grassler Georg
"""

import numpy as np
from scipy.interpolate import CubicSpline

class FluxLableConverter:
    """
    This is a class for the conversion between the poloidal and toroidal flux lable s_pol and s_tor of a given flux surface.
    The conversion is based on the safety factor profile q as q = d{toroidal_flux}/d{polodial_flux}.
    Therefore {toroidal_flux} = int_-inf^{polodial_flux}q({flux})d{flux} - , with {s} = {flux}/{flux_max}
    gives {toroidal_flux} = {flux_max}*int_-inf^{s_pol}q(s)ds
    and {s_tor} = {todroidal_flux}/{toroidal_flux_max} = int_-inf^{s_pol}q(s)ds / int_-inf^1q(s)ds
    """

    def __init__(self,q_profile:np.ndarray):
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

    def spol2stor(self,spol):
        """
        Converts the poloidal flux lable s_pol to the toroidal flux lable s_tor.
        """

        stor = (self.interp_psitor(spol)-self.psitor_min)/(self.psitor_max - self.psitor_min)

        return stor

    def stor2spol(self,stor):
        """
        Converts the toroidal flux lable s_tor to the poloidal flux lable s_pol.
        """

        spol = self.interp_spol(stor)

        return spol

################################################################################
################################################################################
################################################################################

if __name__=='__main__':
    print('This file contains a class for the conversion between the poloidal and') 
    print('toroidal flux lable s_pol and s_tor of a given flux surface.')
    print('Running test conversion for FluxLableConverter.spol2stor:')

    ### TEST conversion

    test_q_profile = np.array([0.107328883E+01, 0.107720493E+01, 0.108176034E+01, 0.108662293E+01, 0.109180789E+01,
    0.109715777E+01, 0.110237995E+01, 0.110759326E+01, 0.111264174E+01, 0.111749784E+01,
    0.112223282E+01, 0.112679797E+01, 0.113119024E+01, 0.113545534E+01, 0.113960896E+01,
    0.114366003E+01, 0.114762476E+01, 0.115151626E+01, 0.115534554E+01, 0.115912742E+01,
    0.116286988E+01, 0.116657936E+01, 0.117026008E+01, 0.117391590E+01, 0.117755199E+01,
    0.118117295E+01, 0.118478192E+01, 0.118838122E+01, 0.119197425E+01, 0.119556509E+01,
    0.119915627E+01, 0.120274966E+01, 0.120634518E+01, 0.120994179E+01, 0.121353789E+01,
    0.121713202E+01, 0.122072380E+01, 0.122431461E+01, 0.122790918E+01, 0.123151371E+01,
    0.123513878E+01, 0.123880044E+01, 0.124251869E+01, 0.124631861E+01, 0.125022176E+01,
    0.125424822E+01, 0.125841538E+01, 0.126274398E+01, 0.126724630E+01, 0.127193166E+01,
    0.127680429E+01, 0.128186393E+01, 0.128710787E+01, 0.129253127E+01, 0.129812616E+01,
    0.130388491E+01, 0.130980104E+01, 0.131586896E+01, 0.132208448E+01, 0.132844437E+01,
    0.133494602E+01, 0.134158692E+01, 0.134836506E+01, 0.135527848E+01, 0.136232650E+01,
    0.136950945E+01, 0.137682756E+01, 0.138428134E+01, 0.139187095E+01, 0.139959587E+01,
    0.140745565E+01, 0.141544973E+01, 0.142357719E+01, 0.143183796E+01, 0.144023278E+01,
    0.144876269E+01, 0.145742867E+01, 0.146623149E+01, 0.147517216E+01, 0.148425197E+01,
    0.149347239E+01, 0.150283494E+01, 0.151234166E+01, 0.152199499E+01, 0.153179756E+01,
    0.154175224E+01, 0.155186207E+01, 0.156213027E+01, 0.157256021E+01, 0.158315538E+01,
    0.159391939E+01, 0.160485603E+01, 0.161596937E+01, 0.162726368E+01, 0.163874355E+01,
    0.165041380E+01, 0.166227938E+01, 0.167434543E+01, 0.168661718E+01, 0.169909986E+01,
    0.171179916E+01, 0.172472097E+01, 0.173787132E+01, 0.175125670E+01, 0.176488396E+01,
    0.177876015E+01, 0.179289220E+01, 0.180728716E+01, 0.182195244E+01, 0.183689591E+01,
    0.185212573E+01, 0.186765021E+01, 0.188347837E+01, 0.189961968E+01, 0.191608329E+01,
    0.193287846E+01, 0.195001494E+01, 0.196750347E+01, 0.198535512E+01, 0.200358108E+01,
    0.202219342E+01, 0.204120481E+01, 0.206062711E+01, 0.208047288E+01, 0.210075619E+01,
    0.212149212E+01, 0.214269631E+01, 0.216438465E+01, 0.218657330E+01, 0.220927901E+01,
    0.223251944E+01, 0.225631341E+01, 0.228068235E+01, 0.230564851E+01, 0.233123430E+01,
    0.235746252E+01, 0.238435665E+01, 0.241194300E+01, 0.244025012E+01, 0.246931034E+01,
    0.249915631E+01, 0.252982035E+01, 0.256133725E+01, 0.259374468E+01, 0.262708650E+01,
    0.266140964E+01, 0.269676294E+01, 0.273319708E+01, 0.277076666E+01, 0.280953348E+01,
    0.284956843E+01, 0.289094597E+01, 0.293374411E+01, 0.297805206E+01, 0.302396718E+01,
    0.307160267E+01, 0.312108021E+01, 0.317253556E+01, 0.322611946E+01, 0.328201171E+01,
    0.334041338E+01, 0.340155473E+01, 0.346569273E+01, 0.353311337E+01, 0.360412353E+01,
    0.367903338E+01, 0.375799260E+01, 0.384097592E+01, 0.392733770E+01, 0.401591803E+01,
    0.410418524E+01, 0.418886581E+01, 0.426572874E+01, 0.433064863E+01, 0.438171292E+01,
    0.441836509E+01, 0.444296595E+01, 0.445858607E+01, 0.446785744E+01, 0.447329187E+01,
    0.447746807E+01])

    converter = FluxLableConverter(test_q_profile)

    test_spol = np.linspace(0.0, 1.0, len(test_q_profile))
    result_stor = converter.spol2stor(test_spol)
    result_spol = converter.stor2spol(result_stor)

    from numpy.testing import assert_allclose
    assert_allclose(result_spol, test_spol)
    print("Interpolation of spol over stor is consistent")
    print("s_pol = ", test_spol[:5],'...',test_spol[-5:])
    print("s_tor = ", result_stor[:5],'...',result_stor[-5:])