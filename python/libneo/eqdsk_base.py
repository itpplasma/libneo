import numpy as np

def readblock(f, k):
    if k > 0 and (k % 5) == 0:
        dat = f.read(17).replace('\n', '')
    else:
        dat = f.read(16)
    return dat

def read_eqdsk(filename):
    eqdata = {}
    with open(filename) as f:
        # First line
        eqdata['header'] = f.read(48)
        line = f.readline().split()
        eqdata['nrgr'] = int(line[1])
        eqdata['nzgr'] = int(line[2])

        # Second line
        eqdata['rboxlength'] = float(f.read(16))
        eqdata['zboxlength'] = float(f.read(16))
        eqdata['R0'] = float(f.read(16))
        eqdata['rboxleft'] = float(f.read(16))
        eqdata['zboxmid'] = float(f.read(16))
        eqdata['R'] = ( eqdata['rboxleft']
                        + np.linspace(0, eqdata['rboxlength'], eqdata['nrgr']) )
        eqdata['Z'] = ( eqdata['zboxmid'] - 0.5*eqdata['zboxlength']
                        + np.linspace(0, eqdata['zboxlength'], eqdata['nzgr']) )
        f.readline()

        # Third line
        # Read magnetic axis
        eqdata['Rpsi0'] = float(f.read(16))  # assumed to be the same as RMAXIS
        eqdata['Zpsi0'] = float(f.read(16))  # assumed to be the same as ZMAXIS
        # Read psiaxis and psiedge
        eqdata['PsiaxisVs'] = float(f.read(16))
        eqdata['PsiedgeVs'] = float(f.read(16))
        # Read B0 and current
        eqdata['Btor_at_R0'] = float(f.read(16))
        f.readline()

        # Fourth line
        eqdata['Ip'] = float(f.read(16))
        f.readline()

        # Fifth line
        f.readline()

        # From sixth line:

        eqdata['fprof'] = np.empty(eqdata['nrgr'])
        for k in range(eqdata['nrgr']):
            eqdata['fprof'][k] = readblock(f, k)
        f.readline()
        eqdata['PsiVs'] = eqdata['PsiVs'].reshape(eqdata['nzgr'], eqdata['nrgr'])

        eqdata['ptotprof'] = np.empty(eqdata['nrgr'])
        for k in range(eqdata['nrgr']):
            eqdata['ptotprof'][k] = readblock(f, k)
        f.readline()

        eqdata['fdfdpsiprof'] = np.empty(eqdata['nrgr'])
        for k in range(eqdata['nrgr']):
            eqdata['fdfdpsiprof'][k] = readblock(f, k)
        f.readline()

        eqdata['dpressdpsiprof'] = np.empty(eqdata['nrgr'])
        for k in range(eqdata['nrgr']):
            eqdata['dpressdpsiprof'][k] = readblock(f, k)
        f.readline()

        eqdata['PsiVs'] = np.empty(eqdata['nrgr']*eqdata['nzgr'])
        for k in range((eqdata['nrgr'])*(eqdata['nzgr'])):
            eqdata['PsiVs'][k] = readblock(f, k)
        f.readline()
        eqdata['PsiVs'] = eqdata['PsiVs'].reshape(eqdata['nzgr'], eqdata['nrgr'])

        eqdata['qprof'] = np.empty(eqdata['nrgr'])
        eqdata['s_pol'] = np.empty(eqdata['nrgr'])
        for k in range(eqdata['nrgr']):
            eqdata['qprof'][k] = readblock(f, k)
            eqdata['s_pol'][k] = float(k)/float(eqdata['nrgr']-1)
        f.readline()

        line = f.readline().split()
        eqdata['npbound'] = int(line[0])
        eqdata['nplimiter'] = int(line[1])

        eqdata['Lcfs'] = np.empty(2*eqdata['npbound'])
        for k in range(2*eqdata['npbound']):
            eqdata['Lcfs'][k] = readblock(f, k)
        f.readline()
        eqdata['Lcfs'] = eqdata['Lcfs'].reshape(-1, 2)

        eqdata['Limiter'] = np.empty(2*eqdata['nplimiter'])
        for k in range(2*eqdata['nplimiter']):
            eqdata['Limiter'][k] = readblock(f, k)
        f.readline()
        eqdata['Limiter'] = eqdata['Limiter'].reshape(-1, 2)

        return eqdata
