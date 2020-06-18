#%% Init
from numpy import linspace, meshgrid, zeros_like, arctan2, ravel, sqrt, float64, frombuffer, array, cos, sin
import matplotlib.pyplot as plt
from libmc_fffi import libmc_efit, field_eq, parmot_mod
from scipy.interpolate import SmoothBivariateSpline
from scipy.integrate import solve_ivp
from IPython import get_ipython

#r = 0.5
#p = 0.5
#z = 0.5
#
#Br = libmc_efit.new('double', 0.0)
#Bp = libmc_efit.new('double', 0.0)
#Bz = libmc_efit.new('double', 0.0)
#dBrdR = libmc_efit.new('double', 0.0)
#dBrdp = libmc_efit.new('double', 0.0)
#dBrdZ = libmc_efit.new('double', 0.0)
#dBpdR = libmc_efit.new('double', 0.0)
#dBpdp = libmc_efit.new('double', 0.0)
#dBpdZ = libmc_efit.new('double', 0.0)
#dBzdR = libmc_efit.new('double', 0.0)
#dBzdp = libmc_efit.new('double', 0.0)
#dBzdZ = libmc_efit.new('double', 0.0)
#
##%% Testing
#
#psi_pol = field_eq.psif - field_eq.psib
#
#print(f'Br before: {Br[0]}')
#print(f'psi_pol before: {psi_pol}')
#
#libmc_efit.field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
#
#psi_pol = field_eq.psif - field_eq.psib
#
#print(f'Br after: {Br[0]}')
#print(f'psi_pol after: {psi_pol}')
#
#r0 = linspace(80,250,30)
#z0 = linspace(-100,100,35)
#RR, ZZ = meshgrid(r0,z0)
#PSI = zeros_like(RR)
#APHI = zeros_like(RR)
#THETA = zeros_like(RR)
#raxis = 100
#for k in range(len(r0)):
#    r = r0[k]
#    for j in range(len(z0)):
#        z = z0[j]
#        libmc_efit.field(r,0.0,z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
#        psi_pol = field_eq.psif - field_eq.psib
#        PSI[j,k] = -psi_pol
#        APHI[j,k] = psi_pol/r
#        THETA[j,k] = arctan2(z,r-raxis)
#
#fig = plt.figure()
#plt.pcolor(RR, ZZ, PSI, cmap='plasma')
#plt.title('Psi'); plt.xlabel('r / cm'); plt.ylabel('z / cm')
#plt.colorbar()
#plt.show()
#
#plt.pcolor(RR,ZZ,APHI, cmap='plasma')
#plt.colorbar()
#plt.show()
#spl = SmoothBivariateSpline(ravel(RR),ravel(ZZ),ravel(APHI))
#
#""" ############################## Constants ############################## """
#
c = 2.9979e10           # cm/s
qe = 4.8032e-10   # franklin ( = 3.336e-10C)
e_mass = 9.1094e-28     # g
p_mass = 1.6726e-24     # g 
ev = 1.6022e-12         # erg ( = 1e-7J)
am = 2                  # Atomic mass 2 of deuterium ions
Zb = 1                    # Atomic charge 1 of deuterium ions
tempi1 = 0.17549561306e4  # ion temperature
v0 = sqrt(2.0*tempi1*ev/(am*p_mass))        # Reference (thermal) velocity

#""" ########################### Non-GC Orbits ############################ """

## Initial conditions
#R0 = 170.  # cm
#phi0 = 0.0
#Z0 = -20.   # cm
#
#vR0 = v0/sqrt(2)
#vph0 = 0#v0/sqrt(3)
#vZ0 = v0/sqrt(2)
#
#libmc_efit.field(R0, phi0, Z0, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
#psi_pol = field_eq.psif - field_eq.psib
#Aphi = -psi_pol/R0
#
#m = am*p_mass
#pphi = m*vph0 + qe/c*Aphi
#
#def ydot(t,y):
#    """ y = [R, p, Z, vR, vZ] """
#    
#    libmc_efit.field(y[0], y[1], y[2], Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
#    psi_pol = field_eq.psif - field_eq.psib
#    Aphi = -psi_pol/y[0]
#    dAphidr = spl.__call__(y[0], y[2], dx=1)
#    dAphidz = spl.__call__(y[0], y[2], dy=1)
#    
#    Rdot = y[3]
#    phidot = (pphi - qe/c*Aphi)/m
#    Zdot = y[4]
#    vRdot = 1/m*(-qe/c*y[4]*Bp[0] + (pphi - qe/c*Aphi)/m*(pphi/y[0] + qe/c*dAphidr))
#    vZdot = 1/m*(qe/c*y[3]*Bp[0] + (pphi - qe/c*Aphi)/m*qe/c*dAphidz)
#   
#    return [Rdot, phidot, Zdot, vRdot, vZdot]
#
##y0 = array([R0, phi0, Z0, vR0, vZ0])
##
##integ = solve_ivp(ydot, (0,1.3e-4), y0, max_step=1e-5)
##sol = integ.y
#
##get_ipython().run_line_magic('matplotlib', 'qt')
##plt.plot(sol[0,::2], sol[2,::2])
##plt.plot(sol[0,0],sol[2,0],'rx')
##plt.plot(sol[0,-1],sol[2,-1],'ro')
##plt.title('Non-Averaged Orbits')
##plt.ylabel('z / cm')
##plt.xlabel('r / cm')
##plt.show()


""" ############################try############################ """

x = libmc_efit._ffi.new('double[3]')
bmod = libmc_efit.new('double', 0.0)
sqrtg = libmc_efit.new('double', 0.0)
bder = libmc_efit._ffi.new('double[3]')
hcovar = libmc_efit._ffi.new('double[3]')
hctrvr = libmc_efit._ffi.new('double[3]')
hcurl = libmc_efit._ffi.new('double[3]')

x[0] = 1e-8      # R0
x[1] = 0.0       # phi0
x[2] = 1e-8      # Z0

libmc_efit.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)


# Reference Larmor radius of thermal particles
bmod_ref = 1e4            # 1 Tesla in Gauss
bmod00 = 1.0              # 1 Tesla in Tesla
rlarm = v0*am*p_mass*c/(Zb*qe*bmod_ref)  # Larmor radius in bmod_ref
parmot_mod.ro0 = rlarm*bmod00                  # Rescaled Larmor radius

# Inverse relativistic temperature, set >=1e5 to neglect relativistic effects
parmot_mod.rmu = 1e5

z = libmc_efit._ffi.new('double[5]')
z[0] = 170.    # R0
z[1] = 0.0     # phi0
z[2] = 20.     # Z0
z[3] = 1.      # normalized velocity module v / v_0
z[4] = 0.      # pitch v_\parallel / v:

x[0:3] = z[0:3]
libmc_efit.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
z0 = frombuffer(libmc_efit._ffi.buffer(z), dtype=float64)

tau = libmc_efit.new('double', 0.0)
vz = libmc_efit._ffi.new('double[5]')  


def velo(t,y):
    z[0:5] = y[:]
    libmc_efit.velo(tau, z, vz)
    return frombuffer(libmc_efit._ffi.buffer(vz), dtype=float64)

dtau = 1.
nt = 10000
times = linspace(0, nt*dtau, nt)
sol = solve_ivp(velo, (times[0], times[-1]), z0, method='LSODA', max_step = 10)
zs = sol.y


