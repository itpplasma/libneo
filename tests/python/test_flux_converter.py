"""
Test for flux_label_converter.py
"""
# %% standard modules
import numpy as np
from numpy.testing import assert_allclose
import pytest

# module to test
from libneo import FluxConverter

# testing data
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

test_axis_polflux = -133.14471161454674
test_edge_polflux = 0.0

def test_FluxConverter():

    test_polflux_profile = np.linspace(test_axis_polflux, test_edge_polflux, len(test_q_profile))
    converter = FluxConverter(test_q_profile, test_polflux_profile)

    result_torflux_profile = converter.polflux2torflux(test_polflux_profile)
    result_polflux_profile = converter.torflux2polflux(result_torflux_profile)

    assert_allclose(result_polflux_profile, test_polflux_profile, atol=1e-15) # choice of atol due to used control data precision
    print("Interpolation of polflux over torflux is selfconsistent")

    control_torflux_profile = np.array([  0.        ,   0.79529968,   1.59375918,   2.39571   ,
    3.2013764 ,   4.01095914,   4.82445535,   5.64180945,
    6.46296881,   7.28778702,   8.11615171,   8.94795915,
    9.78307726,  10.62139507,  11.46282612,  12.30729084,
    13.15471918,  14.00505247,  14.85824041,  15.71424257,
    16.57302718,  17.43456753,  18.2988409 ,  19.16582737,
    20.03551042,  20.90787722,  21.78291783,  22.66062423,
    23.54099043,  24.42401336,  25.30969239,  26.19802855,
    27.08902355,  27.98267862,  28.87899399,  29.7779687 ,
    30.67960109,  31.58388973,  32.49083552,  33.4004434 ,
    34.31272427,  35.22769903,  36.14540177,  37.06588361,
    37.98921312,  38.91547426,  39.84476456,  40.77719592,
    41.71289275,  42.65198716,  43.59461646,  44.54091928,
    45.49103306,  46.44509253,  47.40322759,  48.36556218,
    49.33221506,  50.3033005 ,  51.27892915,  52.25920894,
    53.24424569,  54.23414329,  55.22900401,  56.2289286 ,
    57.23401677,  58.2443682 ,  59.26008276,  60.28126063,
    61.30800226,  62.34040795,  63.37857761,  64.42261081,
    65.47260655,  66.52866339,  67.59088018,  68.65935646,
    69.73419252,  70.8154893 ,  71.90334838,  72.99787218,
    74.09916417,  75.2073289 ,  76.32247221,  77.44470164,
    78.5741266 ,  79.71085853,  80.85501105,  82.00670011,
    83.16604409,  84.33316391,  85.50818312,  86.691228  ,
    87.88242772,  89.08191457,  90.28982407,  91.50629526,
    92.73147078,  93.965497  ,  95.20852416,  96.46070631,
    97.72220158,  98.99317235, 100.27378543, 101.56421225,
    102.86462916, 104.17521769, 105.49616455, 106.82766157,
    108.16990595, 109.52310052, 110.88745406, 112.26318145,
    113.65050392, 115.04964961, 116.46085354, 117.88435754,
    119.32041041, 120.76926853, 122.23119641, 123.70646675,
    125.19536084, 126.69816916, 128.21519132, 129.7467358 ,
    131.29312098, 132.85467603, 134.43174151, 136.02466965,
    137.63382453, 139.25958237, 140.90233211, 142.56247611,
    144.24043163, 145.93663219, 147.6515278 , 149.38558524,
    151.13928825, 152.913139  , 154.70765986, 156.52339579,
    158.36091598, 160.22081353, 162.10370626, 164.01023867,
    165.94108538, 167.89695477, 169.87859069, 171.88677372,
    173.92232321, 175.98610133, 178.07901947, 180.2020431 ,
    182.35619366, 184.54255438, 186.76227728, 189.01659229,
    191.3068163 , 193.63436124, 196.00074442, 198.4076052 ,
    200.85672394, 203.35004022, 205.8896738 , 208.47794577,
    211.11739546, 213.8107988 , 216.56110283, 219.37132098,
    222.24422111, 225.18193541, 228.18523867, 231.25273699,
    234.38026304, 237.56041553, 240.78354519, 244.03906112,
    247.31703819, 250.60971409, 253.91144158, 257.21845125,
    260.52885531])

    assert_allclose(result_torflux_profile, control_torflux_profile)
    print("Interpolation of torflux over polflux is consistent with control data")

    print('----------------------------------------------------------------')
    print('Result of conversion:')
    print("polflux = ", test_polflux_profile[:5],'...',test_polflux_profile[-5:])
    print("torflux = ", result_torflux_profile[:5],'...',result_torflux_profile[-5:])

def can_import_FluxLabelConverter():
    try:
        from libneo import FluxLabelConverter
        return True
    except ImportError:
        return False

@pytest.mark.skipif(not can_import_FluxLabelConverter(), reason="FluxLabelConverter not available")
def test_compare_FluxConverter_to_FluxLabelConverter():

    from libneo import FluxLabelConverter

    # standard calculation of torflux profile with FluxConverter
    test_polflux_profile = np.linspace(test_axis_polflux, test_edge_polflux, len(test_q_profile))
    converter = FluxConverter(test_q_profile, test_polflux_profile)
    result_torflux_profile = converter.polflux2torflux(test_polflux_profile)

    # alternative calculation of torflux profile with FluxLabelConverter
    label_converter = FluxLabelConverter(test_q_profile)
    s_pol = np.linspace(0.0, 1.0, len(test_q_profile))
    s_tor = label_converter.spol2stor(s_pol)
    alternative_torflux_profile = s_tor * label_converter.torflux_max_over_delta_polflux * (test_polflux_profile[-1]-test_polflux_profile[0])

    assert_allclose(alternative_torflux_profile, result_torflux_profile)
    print("Interpolation of torflux over polflux is consistent with alternative FluxLabelConverter result")