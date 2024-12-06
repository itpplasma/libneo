import numpy as np
import matplotlib.pyplot as plt

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


def plot_poloidal_flux(eqdsks: list, labels: list, n_levels: int=10):
    eqdsks = make_to_list_if_not(eqdsks)
    labels = make_to_list_if_not(labels)
    n_eqdsk = len(eqdsks)
    for idx in range(n_eqdsk):
        eqdsk = eqdsks[idx]
        label = labels[idx]
        plt.figure()
        levels = np.linspace(eqdsk["PsiaxisVs"], eqdsk["PsiedgeVs"], n_levels)
        if (levels[1]-levels[0]) < 0:
            levels = levels[-1::-1]
        plt.contourf(eqdsk["R"], eqdsk["Z"], eqdsk["PsiVs"], levels=levels)
        plt.colorbar()
        plt.axis("equal")
        plt.xlabel(r"$R \mathrm{ [m]}$")
        plt.ylabel(r"$Z \mathrm{ [m]}$")
        plt.title(label)
        plt.grid(True)
    plt.show()


def plot_fpol(eqdsks: list, labels: list, markers: list=[]):
    eqdsks = make_to_list_if_not(eqdsks)
    labels = make_to_list_if_not(labels)
    n_eqdsk = len(eqdsks)
    _, ax = plt.subplots(nrows=max(n_eqdsk,2), ncols=1)
    if len(markers)==0:
        markers = ["-b"] * n_eqdsk
    for idx in range(n_eqdsk):
        eqdsk = eqdsks[idx]
        label = labels[idx]
        marker = markers[idx]
        sqrtspol= np.sqrt(eqdsk["s_pol"])
        ax[idx].plot(sqrtspol, eqdsk["fprof"], marker, label=label + " / " r"$B_\mathrm{toroidal}$ at axis=" + str(eqdsk["Btor_at_R0"]))
        ax[idx].set_ylabel(r"$f_\mathrm{pol}$")
        ax[idx].legend()
    ax[-1].set_xlabel(r"$\rho_\mathrm{pol}$")
    plt.show()


def plot_fluxsurfaces(eqdsks: list, labels: list, n_surf: int=11):
    eqdsks = make_to_list_if_not(eqdsks)
    labels = make_to_list_if_not(labels)
    n_eqdsk = len(eqdsks)
    for idx in range(n_eqdsk):
        eqdsk = eqdsks[idx]
        label = labels[idx]
        normalized_flux = ((eqdsk["PsiVs"] - eqdsk["PsiaxisVs"]) / 
                           (eqdsk["PsiedgeVs"] - eqdsk["PsiaxisVs"]))
        plt.figure()
        levels = np.linspace(0, 1, n_surf)
        fluxsurfaces = plt.contour(eqdsk["R"], eqdsk["Z"], normalized_flux, levels=levels)
        plt.clabel(fluxsurfaces, inline=True, fontsize=8)
        plt.axis("equal")
        plt.xlabel(r"$R \mathrm{ [m]}$")
        plt.ylabel(r"$Z \mathrm{ [m]}$")
        plt.title(label)
        plt.grid(True)
    plt.show()


def make_to_list_if_not(obj):
    if not isinstance(obj, list):
        obj = [obj]
    return obj
