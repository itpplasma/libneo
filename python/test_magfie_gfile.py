#%% Init
from numpy import linspace, sqrt, meshgrid, cos, sin, zeros_like, arctan2
import matplotlib.pyplot as plt
from libmc_fffi import libmc_efit, field_eq

r = libmc_efit.new('double', 50.0)
p = libmc_efit.new('double', 0.0)
z = libmc_efit.new('double', 0.0)

Br = libmc_efit.new('double', 0.0)
Bp = libmc_efit.new('double', 0.0)
Bz = libmc_efit.new('double', 0.0)
dBrdR = libmc_efit.new('double', 0.0)
dBrdp = libmc_efit.new('double', 0.0)
dBrdZ = libmc_efit.new('double', 0.0)
dBpdR = libmc_efit.new('double', 0.0)
dBpdp = libmc_efit.new('double', 0.0)
dBpdZ = libmc_efit.new('double', 0.0)
dBzdR = libmc_efit.new('double', 0.0)
dBzdp = libmc_efit.new('double', 0.0)
dBzdZ = libmc_efit.new('double', 0.0)

#%% Testing

psi_pol = field_eq.psif - field_eq.psib

print(f'Br before: {Br[0]}')
print(f'psi_pol before: {psi_pol}')

libmc_efit.field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

psi_pol = field_eq.psif - field_eq.psib

print(f'Br after: {Br[0]}')
print(f'psi_pol after: {psi_pol}')

# check EFIT file with https://github.com/itpplasma/libneo/tree/master/matlab/EFIT

#%% spline psi

# grad A_(phi) = grad A_phi/R = grad (-psi/R)

# RR, ZZ = meshgrid(...)

r0 = linspace(70, 250, 60)
z0 = linspace(-120, 120, 60)
RR, ZZ = meshgrid(r0,z0)
PSI = zeros_like(RR)
PHI = zeros_like(RR)
for k in range(len(r0)):
    r[0] = r0[k]
    for j in range(len(z0)):
        z[0] = z0[j]
        libmc_efit.field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
        psi_pol = field_eq.psif - field_eq.psib
        PSI[j,k] = -psi_pol
        
        PHI[j,k] = arctan2(z[0], r[0])
#fig, axs = plt.figure()
        
fig = plt.figure()
plt.contour(RR, ZZ, PSI, cmap='plasma', levels=50)
#fig, ax = plt.subplots()
#CS = ax.contour(RR, ZZ, PSI)
#ax.clabel(CS, inline=1, fontsize=10)
#ax.set_title('Simplest default with labels')
plt.title('Magnetic flux'); plt.xlabel('r / cm'); plt.ylabel('z / cm')
plt.colorbar()
plt.show()
# for ...:
#    PSI[k,l] = -psi/R

# psifun = spline(RR, ZZ, PSI)
# e.g. scipy.interpolate.SmoothBivariateSpline

# plot contours of psi
