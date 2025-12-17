import numpy as np
from scipy.interpolate import PchipInterpolator, LinearNDInterpolator, interp1d
from scipy.integrate import cumulative_trapezoid, quad
from datetime import datetime
from netCDF4 import Dataset
from libneo import eqdsk_base

def read_vmec_data(filename):
    dataset = Dataset(filename, mode='r')
    data = {}

    data['global_attributes'] = {attr: getattr(dataset, attr) for attr in dataset.ncattrs()}

    data['dimensions'] = {dim: dataset.dimensions[dim].size for dim in dataset.dimensions}

    variables = {}
    for var_name in dataset.variables:
        var = dataset.variables[var_name]
        var_data = var[:]
        var_attrs = {attr: getattr(var, attr) for attr in var.ncattrs()}
        variables[var_name] = {
            'data': var_data,
            'attributes': var_attrs
        }
    data['variables'] = variables

    dataset.close()
    return data

def vmec_to_eqdsk(vmec_data):
    
    ns = int(vmec_data['variables']['ns']['data'])
    nx = ns
    nz = ns
    ntheta = 360
    theta = np.linspace(0, 2*np.pi, ntheta)
    zeta = 0.0  # For axisymmetric output

    rmnc = vmec_data['variables']['rmnc']['data'].T
    zmns = vmec_data['variables']['zmns']['data'].T
    xm    = vmec_data['variables']['xm']['data']
    xn    = vmec_data['variables']['xn']['data']

    r = cfunct(theta, zeta, rmnc, xm, xn)
    z = sfunct(theta, zeta, zmns, xm, xn)

    iasym = bool(vmec_data['variables'].get('lasym__logical__', {}).get('data', 0))
    if iasym:
        rmns = vmec_data['variables']['rmns']['data'].T  # shape (nmode, ns)
        zmnc = vmec_data['variables']['zmnc']['data'].T  # shape (nmode, ns)
        r += sfunct(theta, zeta, rmns, xm, xn)
        z += cfunct(theta, zeta, zmnc, xm, xn)

    phi = vmec_data['variables']['phi']['data']
    phi=phi/(2*np.pi) #divide by 2pi to get flux in Weber/rad as in eqdsk
    phin = phi / phi[-1]


    """iotaf = h2f(vmec_data['variables']['iotas']['data'], ns)
    chipf = vmec_data['variables']['q_factor']['data'] * vmec_data['variables']['phipf']['data']
    chif = np.cumsum(chipf)
    #print(chif.shape)
    #chif = rho_pol
    chirz = np.tile(chif[:, None], (1, ntheta))"""
    
    # new integration for toroidal flux to poloidal flux
    print(phi)
    iota_interp=interp1d(phi, vmec_data['variables']['iotas']['data'], kind='linear', fill_value="extrapolate")
    chif= np.zeros(len(phi))
    for i in range(len(phi)):
        chif[i], _ = quad(iota_interp, phi[0], phi[i])
    chirz = np.tile(chif[:, None], (1, ntheta))

    points = np.column_stack((r.flatten(), z.flatten()))
    Fchi = LinearNDInterpolator(points, chirz.flatten())

    rminsurf = vmec_data['variables']['rmin_surf']['data']
    rmaxsurf = vmec_data['variables']['rmax_surf']['data']
    zmaxsurf = np.max(z[-1, :]) * 1.1
    zminsurf = np.min(z[-1, :]) * 1.1

    rboxlength = (rmaxsurf - rminsurf) * 1.1
    zboxlength = (zmaxsurf - zminsurf) * 1.1
    R0 = 0.5 * (rmaxsurf + rminsurf)
    rboxleft = R0 - 0.5 * rboxlength
    zboxmid = 0.5 * (zmaxsurf + zminsurf)

    r1 = rboxleft
    r2 = rboxleft + rboxlength
    z1 = zboxmid - 0.5 * zboxlength
    z2 = zboxmid + 0.5 * zboxlength

    xgrid, zgrid = np.meshgrid(np.linspace(r1, r2, nx),
                               np.linspace(z1, z2, nz))
    psixz = Fchi(xgrid, zgrid)
    psixz = np.where(np.isnan(psixz), 0.0, psixz)

    eqdata = {}
    eqdata['header'] = f"{'  VMEC':<10}{datetime.now().strftime('%m/%d/%Y'):>10}{'  #000001':>10}{'  0000ms':>10}       "

    eqdata['nrgr'] = nx
    eqdata['nzgr'] = nz

    eqdata['rboxlength'] = rboxlength
    eqdata['zboxlength'] =  zboxlength
    eqdata['R0'] = R0
    eqdata['rboxleft'] = rboxleft
    eqdata['zboxmid'] = zboxmid

    eqdata['Rpsi0'] = r[0, 0]
    eqdata['Zpsi0'] = z[0, 0]
    eqdata['PsiaxisVs'] = chif[0]
    eqdata['PsiedgeVs'] = chif[-1]
    b0 = vmec_data['variables']['b0']['data']
    eqdata['Btor_at_R0'] = b0

    Itor = vmec_data['variables']['ctor']['data']
    eqdata['Ip'] = Itor

    xnew = np.linspace(0, 1, nx)
    eqdata['fprof'] = PchipInterpolator(phin, phi)(xnew)
    presf = vmec_data['variables']['presf']['data']
    eqdata['ptotprof'] = PchipInterpolator(phin, presf)(xnew)
    phi_grad = np.gradient(phi)
    eqdata['fdfdpsiprof'] = PchipInterpolator(phin, phi_grad)(xnew)
    presf_grad = np.gradient(presf)
    eqdata['dpressdpsiprof'] = PchipInterpolator(phin, presf_grad)(xnew)

    eqdata['PsiVs'] = psixz.flatten()

    if 'q_factor' in vmec_data['variables']:
        print("writing q profile")
        q_factor = vmec_data['variables']['q_factor']['data']
        eqdata['qprof'] = PchipInterpolator(np.linspace(0,1,len(q_factor)), q_factor)(np.linspace(0,1,nx))
    else:
        eqdata['qprof'] = np.linspace(1, 2, nx)

    eqdata['npbound'] = ntheta
    eqdata['nplimiter'] = 0
    eqdata['Lcfs'] = np.stack((r[-1,:],z[-1,:]),axis=-1)
    eqdata['Limiter'] = np.empty((0, 2))

    return eqdata

def h2f(var, ns):
    temp = np.zeros(ns)
    temp[0] = 1.5 * var[1] - 0.5 * var[2]
    temp[1:ns-1] = 0.5 * (var[1:ns-1] + var[2:ns])
    temp[ns-1] = 2 * var[ns-1] - var[ns-2]
    return temp


def cfunct(theta, zeta, coeff, xm, xn):
    nmode, ns = coeff.shape
    result = np.zeros((ns, len(theta)))
    for m in range(nmode):
        result += coeff[m, :].reshape(ns, 1) * np.cos(xm[m] * theta - xn[m] * zeta)
    return result

def sfunct(theta, zeta, coeff, xm, xn):
    nmode, ns = coeff.shape
    result = np.zeros((ns, len(theta)))
    for m in range(nmode):
        result += coeff[m, :].reshape(ns, 1) * np.sin(xm[m] * theta - xn[m] * zeta)
    return result
