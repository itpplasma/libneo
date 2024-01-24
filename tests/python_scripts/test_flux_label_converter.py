"""
Test for flux_label_converter.py
"""
# standard modules
import numpy as np

# module to test
from libneo import FluxLabelConverter

def test_FluxLabelConverter():

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

    converter = FluxLabelConverter(test_q_profile)

    test_spol = np.linspace(0.0, 1.0, len(test_q_profile))
    result_stor = converter.spol2stor(test_spol)
    result_spol = converter.stor2spol(result_stor)

    from numpy.testing import assert_allclose
    assert_allclose(result_spol, test_spol)
    print("Interpolation of spol over stor is selfconsistent")

    control_stor = np.array([0.        , 0.00305264, 0.0061174 , 0.00919556, 0.01228799,
       0.01539545, 0.01851793, 0.02165522, 0.02480711, 0.02797305,
       0.0311526 , 0.03434537, 0.03755084, 0.04076859, 0.0439983 ,
       0.04723965, 0.05049237, 0.05375624, 0.05703107, 0.06031671,
       0.06361302, 0.06691991, 0.07023729, 0.07356508, 0.07690323,
       0.08025168, 0.08361038, 0.08697933, 0.09035848, 0.09374782,
       0.09714737, 0.10055711, 0.10397706, 0.10740721, 0.11084758,
       0.11429816, 0.11775894, 0.12122991, 0.12471108, 0.12820247,
       0.13170412, 0.13521611, 0.13873857, 0.1422717 , 0.14581576,
       0.14937107, 0.15293801, 0.156517  , 0.16010853, 0.1637131 ,
       0.16733124, 0.17096348, 0.17461034, 0.17827235, 0.18195001,
       0.18564378, 0.18935413, 0.19308149, 0.19682629, 0.20058895,
       0.20436986, 0.20816943, 0.21198805, 0.21582611, 0.21968398,
       0.22356206, 0.22746073, 0.23138036, 0.23532135, 0.23928408,
       0.24326894, 0.2472763 , 0.25130655, 0.25536006, 0.25943721,
       0.2635384 , 0.26766399, 0.27181438, 0.27598996, 0.28019112,
       0.28441826, 0.28867178, 0.29295209, 0.29725959, 0.30159472,
       0.30595789, 0.31034954, 0.31477012, 0.31922009, 0.3236999 ,
       0.32821003, 0.33275096, 0.3373232 , 0.34192725, 0.34656362,
       0.35123286, 0.35593551, 0.36067213, 0.3654433 , 0.37024961,
       0.37509166, 0.37997009, 0.38488553, 0.38983863, 0.39483008,
       0.39986057, 0.40493083, 0.41004157, 0.41519357, 0.4203876 ,
       0.42562446, 0.43090498, 0.43623001, 0.44160041, 0.4470171 ,
       0.452481  , 0.45799307, 0.46355429, 0.46916568, 0.47482827,
       0.48054317, 0.48631146, 0.49213432, 0.49801292, 0.50394848,
       0.50994227, 0.51599559, 0.5221098 , 0.5282863 , 0.53452652,
       0.54083196, 0.54720417, 0.55364474, 0.56015535, 0.56673771,
       0.57339363, 0.58012495, 0.5869336 , 0.59382159, 0.60079102,
       0.60784406, 0.61498299, 0.62221018, 0.62952811, 0.63693937,
       0.64444668, 0.65205288, 0.65976098, 0.66757413, 0.67549562,
       0.68352897, 0.69167787, 0.69994624, 0.70833825, 0.71685832,
       0.72551116, 0.73430183, 0.74323576, 0.75231876, 0.76155712,
       0.77095769, 0.7805279 , 0.79027589, 0.80021058, 0.8103417 ,
       0.82067991, 0.83123653, 0.84202313, 0.85305031, 0.86432628,
       0.87585399, 0.88762812, 0.89963264, 0.91183917, 0.92421066,
       0.93670646, 0.94928847, 0.9619269 , 0.97460007, 0.98729352,
       1.        ])

    assert_allclose(result_stor, control_stor, atol=1e-6) # choice of atol due to used control data precision
    print("Interpolation of stor over spol is consistent with control data")

    print('----------------------------------------------------------------')
    print('Result of conversion:')
    print("s_pol = ", test_spol[:5],'...',test_spol[-5:])
    print("s_tor = ", result_stor[:5],'...',result_stor[-5:])