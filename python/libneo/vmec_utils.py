import numpy as np
from os.path import splitext
from scipy.interpolate import PchipInterpolator, LinearNDInterpolator, interp1d
from scipy.integrate import quad
from datetime import datetime
from netCDF4 import Dataset
from libneo import eqdsk_base

def eqdsk2vmec(eqdsk_file, vmec_in_file=None):
    data = eqdsk2vmec_gfile(eqdsk_file)

    if vmec_in_file is None:
        vmec_in_file = f'input.{splitext(eqdsk_file)[0]}'

    # Initialize VMEC input
    vmec_input = default_vmec_input(data)

    eqdsk_file = eqdsk_file.split('.')[0]
    write_vmec_input(vmec_in_file, vmec_input)

def default_vmec_input(data):
    return{
        "delt": 1.0,
        "niter": 20000,
        "tcon0": 1.0,
        "ns_array": [16, 32, 64, 128],
        "ftol_array": [1e-30, 1e-30, 1e-30, 1e-12],
        "niter_array": [1000, 2000, 4000, 20000],
        "lasym": 1,
        "nfp": 1,
        "mpol": len(data['rbc']),
        "ntor": 0,
        "nzeta": 1,
        "nstep": 200,
        "ntheta": 2 * len(data['rbc']) + 6,
        "phiedge": data['phiedge'],
        "lfreeb": 0,
        "mgrid_file": '',
        "extcur": [0, 0, 0],
        "nvacskip": 6,
        "gamma": 0.0,
        "bloat": 1.0,
        "spres_ped": 1.0,
        "pres_scale": 1.0,
        "pmass_type": 'akima_spline',
        "am": data['am'],
        "am_aux_s": data['am_aux_s'],
        "am_aux_f": data['am_aux_f'],
        "pcurr_type": 'akima_spline_ip',
        "curtor": data['curtor'],
        "ncurr": 0,
        "ac": data['ac'],
        "ac_aux_s": data['ac_aux_s'],
        "ac_aux_f": data['ac_aux_f'],
        "piota_type": 'akima_spline',
        "ai": data['ai'],
        "ai_aux_s": data['ai_aux_s'],
        "ai_aux_f": data['ai_aux_f'],
        "raxis_cc": data['raxis'],
        "raxis_cs": 0.0,
        "zaxis_cc": data['zaxis'],
        "zaxis_cs": 0.0,
        "rbc": data['rbc'].T,
        "zbs": data['zbs'].T,
        "rbs": data['rbs'].T,
        "zbc": data['zbc'].T
    }

def eqdsk2vmec_gfile(filename):
    data = eqdsk_base.read_eqdsk(filename)

    # Extract boundary and profiles
    R = data['Lcfs'][:, 0]
    Z = data['Lcfs'][:, 1]
    theta = np.arctan2(Z, R - np.mean(R))
    theta[theta < 0] += 2 * np.pi
    sorted_indices = np.argsort(theta)
    R, Z = R[sorted_indices], Z[sorted_indices]

    # Interpolate flux functions
    pflux = np.linspace(data['PsiaxisVs'], data['PsiedgeVs'], len(data['fprof']))

    press = data['ptotprof']
    qprof = data['qprof']

    qprof_interp = interp1d(pflux, qprof, kind='linear', fill_value="extrapolate")

    torflux = np.zeros(len(pflux))

    for i in range(len(pflux)):
        torflux[i], _ = quad(qprof_interp, pflux[0], pflux[i])   # EQDSK uses flux in Weber/rad

    torflux_vmec = torflux * 2 * np.pi  # VMEC uses flux in Weber
    torflux_norm = torflux/torflux[-1]

    am_aux_s = np.linspace(0, 1, 100)
    am_aux_f = PchipInterpolator(torflux_norm, press)(am_aux_s)
    p = np.polyfit(am_aux_s, am_aux_f, 9)
    p[-1] = am_aux_f[0]
    p = np.concatenate(([-np.sum(p[:10]) + am_aux_f[-1]], p))
    am = np.flip(p)

    jdotb= data['dpressdpsiprof']+data['fdfdpsiprof']
    ac_aux_s = np.linspace(0, 1, 100)
    ac_aux_f = PchipInterpolator(torflux_norm, jdotb)(ac_aux_s)
    p = np.polyfit(ac_aux_s, ac_aux_f, 9)
    p[-1] = ac_aux_f[0]
    p = np.concatenate(([-np.sum(p[:10]) + ac_aux_f[-1]], p))
    ac = np.flip(p)

    ai_aux_s = np.linspace(0, 1, 100)
    ai_aux_f = PchipInterpolator(torflux_norm, 1/data['qprof'])(ai_aux_s)
    p = np.polyfit(ai_aux_s, ai_aux_f, 9)
    p[-1] = ai_aux_f[0]
    p = np.concatenate(([-np.sum(p[:10]) + ac_aux_f[-1]], p))
    ai = np.flip(p)

    data2 = {
        'phiedge': torflux_vmec[-1],
        'curtor': data['Ip'],
        'raxis': data['Rpsi0'],
        'zaxis': data['Zpsi0'],
        'am': am,
        'am_aux_s': am_aux_s,
        'am_aux_f': am_aux_f,
        'ac': ac,
        'ac_aux_s': ac_aux_s,
        'ac_aux_f': ac_aux_f,
        'ai': ai,
        'ai_aux_s': ai_aux_s,
        'ai_aux_f': ai_aux_f
    }

    compute_boundary_fourier(R, Z, data2)

    return data2

def compute_boundary_fourier(R, Z, data, mpol=12, ntor=0):
    data['rbc'] = np.zeros((mpol+1, ntor+1))
    data['zbs'] = np.zeros((mpol+1, ntor+1))
    data['rbs'] = np.zeros((mpol+1, ntor+1))
    data['zbc'] = np.zeros((mpol+1, ntor+1))

    nu = len(R)
    nv = 1

    cosu = np.zeros((nu, mpol + 1))
    cosv = np.zeros((nv, 2 * ntor + 1))
    sinu = np.zeros((nu, mpol + 1))
    sinv = np.zeros((nv, 2 * ntor + 1))

    alu = 2 * np.pi / nu
    for i in range(nu):
        for j in range(mpol + 1):
            m = j
            cosu[i, j] = np.cos(m * i * alu)
            sinu[i, j] = np.sin(m * i * alu)

    alv = 2 * np.pi / nv
    for i in range(nv):
        for j in range(2 * ntor + 1):
            n = j - ntor
            cosv[i, j] = np.cos(n * i * alv)
            sinv[i, j] = np.sin(n * i * alv)

    # Initialize fnuv
    fnuv = np.zeros(mpol + 1)
    fnuv[0] = 1.0 / (nu * nv)
    for i in range(1, mpol + 1):
        fnuv[i] = 2 * fnuv[0]

    for m1 in range(mpol + 1):
        for n1 in range(2 * ntor + 1):
            for i in range(nu):
                for j in range(nv):
                    data['rbc'][m1, n1] += R[i] * (cosv[j, n1] * cosu[i, m1] - sinv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    data['zbs'][m1, n1] += Z[i] * (sinv[j, n1] * cosu[i, m1] + cosv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    data['rbs'][m1, n1] += R[i] * (sinv[j, n1] * cosu[i, m1] + cosv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    data['zbc'][m1, n1] += Z[i] * (cosv[j, n1] * cosu[i, m1] - sinv[j, n1] * sinu[i, m1]) * fnuv[m1]


def write_vmec_input_vec(f, vec):
    for i in range(len(vec)):
        f.write(f"{vec[i]:16.14E} ")
        if i%3==2:
            f.write('\n ')

    f.write('\n')

def write_vmec_input(filename, data):
    """
    Writes a VMEC input structure to a file.
    """
    """    if 'datatype' not in data or data['datatype'] != 'VMEC_input':
        print("Error: Not a valid VMEC input data structure!")
        return -1"""

    try:
        with open(filename, 'w') as f:
            f.write('&INDATA\n')
            f.write('!----- Runtime Parameters -----\n')
            f.write(f'  DELT = {data["delt"]:16.14E}\n')
            f.write(f'  NITER = {data["niter"]}\n')
            f.write(f'  NSTEP = {data["nstep"]}\n')
            f.write(f'  TCON0 = {data["tcon0"]:16.14E}\n')

            ns = len(data['ns_array'])
            f.write('  NS_ARRAY = ' + ' '.join(f'{val:8d}' for val in data['ns_array']) + '\n')

            if 'niter_array' in data:
                f.write('  NITER_ARRAY = ' + ' '.join(f'{val:8d}' for val in data['niter_array']) + '\n')

            f.write('  FTOL_ARRAY = ' + ' '.join(f'{val:.2E}' for val in data['ftol_array']) + '\n')

            f.write('!----- Grid Parameters -----\n')
            f.write(f'  LASYM = {bool(data["lasym"])}\n')
            f.write(f'  NFP = {data["nfp"]}\n')
            f.write(f'  MPOL = {data["mpol"]}\n')
            f.write(f'  NTOR = {data["ntor"]}\n')
            if 'ntheta' in data:
                f.write(f'  NTHETA = {data["ntheta"]}\n')
            if 'nzeta' in data:
                f.write(f'  NZETA = {data["nzeta"]}\n')
            f.write(f'  PHIEDGE = {data["phiedge"]}\n')

            f.write('!----- Free Boundary Parameters -----\n')
            f.write(f'  LFREEB = {bool(data["lfreeb"])}\n')
            if data['lfreeb']:
                f.write(f'  MGRID_FILE = "{data["mgrid_file"]}"\n')
                f.write('  EXTCUR = ' + ' '.join(map(str, data['extcur'])) + '\n')
                f.write(f'  NVACSKIP = {data["nvacskip"]}\n')

            f.write('!----- Pressure Parameters -----\n')
            f.write(f'  GAMMA = {data["gamma"]:16.14E}\n')
            f.write(f'  BLOAT = {data["bloat"]:16.14E}\n')
            f.write(f'  SPRES_PED = {data["spres_ped"]:16.14E}\n')
            f.write(f'  PMASS_TYPE = "{data["pmass_type"]}"\n')

            if 'pres_scale' in data:
                f.write(f'  PRES_SCALE = {data["pres_scale"]:16.14E}\n')
            if 'am' in data:
                f.write('  AM = ' + ' ')
                write_vmec_input_vec(f, data['am'])
            if 'am_aux_s' in data:
                f.write('  AM_AUX_S = ' + ' ')
                write_vmec_input_vec(f, data['am_aux_s'])
            if 'am_aux_f' in data:
                f.write('  AM_AUX_F = ' + ' ')
                write_vmec_input_vec(f, data['am_aux_f'])

            f.write('!----- Current/Iota Parameters -----\n')
            f.write(f'  CURTOR = {data["curtor"]:16.14E}\n')
            f.write(f'  NCURR = {data["ncurr"]}\n')
            f.write(f'  PIOTA_TYPE = "{data["piota_type"]}"\n')

            if 'ai' in data:
                f.write('  AI = ' + ' ')
                write_vmec_input_vec(f, data['ai'])
            if 'ai_aux_s' in data:
                f.write('  AI_AUX_S = ' + ' ')
                write_vmec_input_vec(f, data['ai_aux_s'])
            if 'ai_aux_f' in data:
                f.write('  AI_AUX_F = ' + ' ')
                write_vmec_input_vec(f, data['ai_aux_f'])
            f.write(f'  PCURR_TYPE = "{data["pcurr_type"]}"\n')
            if 'ac' in data:
                f.write('  AC = ' + ' ')
                write_vmec_input_vec(f, data['ac'])
            if 'ac_aux_s' in data:
                f.write('  AC_AUX_S = ' + ' ')
                write_vmec_input_vec(f, data['ac_aux_s'])
            if 'ac_aux_f' in data:
                f.write('  AC_AUX_F = ' + ' ')
                write_vmec_input_vec(f, data['ac_aux_f'])
            f.write('!----- Axis Parameters -----\n')
            if 'raxis' in data:
                f.write(f'  RAXIS = {data["raxis"]:16.14E}\n')
            elif 'raxis_cc' in data:
                f.write(f'  RAXIS_CC = {data["raxis_cc"]:16.14E}\n')
                if data['lasym']:
                    f.write(f'  RAXIS_CS = {data["raxis_cs"]:16.14E}\n')

            if 'zaxis' in data:
                f.write(f'  ZAXIS = {data["zaxis"]:16.14E}\n')
            elif 'zaxis_cc' in data:
                f.write(f'  ZAXIS_CC = {data["zaxis_cc"]:16.14E}\n')
                if data['lasym']:
                    f.write(f'  ZAXIS_CS = {data["zaxis_cs"]:16.14E}\n')

            f.write('!----- Boundary Parameters -----\n')
            for i in range(data['mpol']):
                for j in range(2 * data['ntor'] + 1):
                    if data['rbc'][j, i] != 0.0 or data['zbs'][j, i] != 0.0:
                        f.write(f'  RBC({j-data["ntor"]},{i}) = {data["rbc"][j,i]:.10e}  ZBS({j-data["ntor"]},{i}) = {data["zbs"][j,i]:.10e}\n')
                        if data['lasym']:
                            f.write(f'     RBS({j-data["ntor"]},{i}) = {data["rbs"][j,i]:.10e}  ZBC({j-data["ntor"]},{i}) = {data["zbc"][j,i]:.10e}\n')

            f.write('!----- Created by write_vmec ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' -----\n')
            f.write('/\n')
        return 1
    except Exception as e:
        print(f"Error writing file: {e}")
        return -1

def read_vmec_out_data(filename):
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

    # new integration for toroidal flux to poloidal flux
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
    eqdata['header'] = f"{'VMEC ':<10}{datetime.now().strftime('%m/%d/%Y'):>10}{'         # xxxxx':>10}{'  0000ms':>10}           0  "

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
    print(vmec_data['variables'].keys())
    fpol=vmec_data['variables']['rbtor']['data']
    eqdata['fprof'] = PchipInterpolator(phin, fpol)(xnew)
    presf = vmec_data['variables']['presf']['data']
    eqdata['ptotprof'] = PchipInterpolator(phin, presf)(xnew)
    phi_grad = np.gradient(phi)
    eqdata['fdfdpsiprof'] = PchipInterpolator(phin, phi_grad)(xnew)
    presf_grad = np.gradient(presf)
    eqdata['dpressdpsiprof'] = PchipInterpolator(phin, presf_grad)(xnew)

    eqdata['PsiVs'] = psixz.flatten()

    if 'q_factor' in vmec_data['variables']:
        q_factor = vmec_data['variables']['q_factor']['data']
        eqdata['qprof'] = PchipInterpolator(np.linspace(0,1,len(q_factor)), q_factor)(np.linspace(0,1,nx))
    else:
        eqdata['qprof'] = np.linspace(1, 2, nx)

    eqdata['npbound'] = ntheta
    eqdata['nplimiter'] = 0
    eqdata['Lcfs'] = np.stack((r[-1,:],z[-1,:]),axis=-1)
    eqdata['Limiter'] = np.empty((0, 2))

    return eqdata

def cfunct(theta, zeta, coeff, xm, xn):
    nmode, ns = coeff.shape
    result = np.zeros((ns, len(theta)))
    for m in range(nmode):
        result += coeff[m, :].reshape(ns, 1) * np.cos(xm[m] * theta + xn[m] * zeta)
    return result

def sfunct(theta, zeta, coeff, xm, xn):
    nmode, ns = coeff.shape
    result = np.zeros((ns, len(theta)))
    for m in range(nmode):
        result += coeff[m, :].reshape(ns, 1) * np.sin(xm[m] * theta + xn[m] * zeta)
    return result
