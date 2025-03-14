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

import numpy as np

def write_eqdsk(filename, eqdata):
    def write_field(fobj, value):
        fobj.write(f"{float(value):16.8E}")

    def write_array(fobj, arr):
        arr = np.atleast_1d(arr)
        for i, val in enumerate(arr):
            write_field(fobj, val)
            if (i + 1) % 5 == 0:
                fobj.write("\n")
        if len(arr) % 5 != 0:
            fobj.write("\n")

    with open(filename, 'w') as f:
        f.write(eqdata['header'][:48])

        f.write(f" 0 {eqdata['nrgr']} {eqdata['nzgr']}\n")

        write_field(f, eqdata['rboxlength'])
        write_field(f, eqdata['zboxlength'])
        write_field(f, eqdata['R0'])
        write_field(f, eqdata['rboxleft'])
        write_field(f, eqdata['zboxmid'])
        f.write("\n")

        write_field(f, eqdata['Rpsi0'])
        write_field(f, eqdata['Zpsi0'])
        write_field(f, eqdata['PsiaxisVs'])
        write_field(f, eqdata['PsiedgeVs'])
        write_field(f, eqdata['Btor_at_R0'])
        f.write("\n")

        write_field(f, eqdata['Ip'])
        write_field(f, eqdata['PsiaxisVs'])
        write_field(f, 0.0)
        write_field(f, eqdata['Rpsi0'])
        write_field(f, 0.0)
        f.write("\n")

        write_field(f, eqdata['Zpsi0'])
        write_field(f, 0.0)
        write_field(f, eqdata['PsiedgeVs'])
        write_field(f, 0.0)
        write_field(f, 0.0)
        f.write("\n")

        write_array(f, eqdata['fprof'])
        write_array(f, eqdata['ptotprof'])
        write_array(f, eqdata['fdfdpsiprof'])
        write_array(f, eqdata['dpressdpsiprof'])

        write_array(f, eqdata['PsiVs'].flatten())

        write_array(f, eqdata['qprof'])
        print("writing q profile to file")
        print(eqdata['qprof'])
        f.write(f"{eqdata['npbound']} {eqdata['nplimiter']}\n")

        print("writing Lcfs to file")
        print(eqdata['Lcfs'])
        write_array(f, eqdata['Lcfs'].flatten())
        
        write_array(f, eqdata['Limiter'].flatten())


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
        plt.xlabel(r"$R \; \mathrm{[m]}$")
        plt.ylabel(r"$Z \; \mathrm{[m]}$")
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
        plt.xlabel(r"$R \; \mathrm{[m]}$")
        plt.ylabel(r"$Z \; \mathrm{[m]}$")
        plt.title(label)
        plt.grid(True)
    plt.show()


def make_to_list_if_not(obj):
    if not isinstance(obj, list):
        obj = [obj]
    return obj


def approximate_fluxsurface_area_eqdsk(eqdsk, n_surf: int=11):
    import numpy as np
    import matplotlib.pyplot as plt

    normalized_flux = ((eqdsk["PsiVs"] - eqdsk["PsiaxisVs"]) /
                            (eqdsk["PsiedgeVs"] - eqdsk["PsiaxisVs"]))
    spol = np.linspace(0, 1, n_surf)
    Zmin = np.min(eqdsk["Lcfs"][:,1])
    Zmax = np.max(eqdsk["Lcfs"][:,1])
    I = (eqdsk["Z"] > Zmin) * (eqdsk["Z"] < Zmax)
    fluxcontour = plt.contour(eqdsk["R"], eqdsk["Z"][I], normalized_flux[I], levels=spol)

    def compute_distance(p1, p2):
        return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

    circumference = []
    for flux_contour in fluxcontour.allsegs:
        total_length = 0
        points = np.array(flux_contour[0])
        for idx in range(1, len(points)):
            total_length += compute_distance(points[idx-1], points[idx])
        circumference.append(total_length)
    plt.close()

    area = 2*np.pi*eqdsk["R0"]*np.array(circumference)
    area = area[:-1] # exclude separatrix
    spol = spol[:-1]
    return area, spol
