#%% Init
import numpy as np
import matplotlib.pyplot as plt
from libmc_fffi import libmc_efit, field_eq, parmot_mod
from scipy.interpolate import SmoothBivariateSpline, UnivariateSpline
from scipy.integrate import solve_ivp
from IPython import get_ipython
from random import random
import time

t = time.time()
r = 0.5
p = 0.5
z = 0.5

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
#
#%% Testing

psi_pol = field_eq.psif - field_eq.psib

print(f'Br before: {Br[0]}')
print(f'psi_pol before: {psi_pol}')

libmc_efit.field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

psi_pol = field_eq.psif - field_eq.psib

print(f'Br after: {Br[0]}')
print(f'psi_pol after: {psi_pol}')

r0 = np.linspace(80,250,30)
z0 = np.linspace(-100,100,35)
RR, ZZ = np.meshgrid(r0,z0)
PSI = np.zeros_like(RR)
APHI = np.zeros_like(RR)
THETA = np.zeros_like(RR)
raxis = 170
for k in range(len(r0)):
    r = r0[k]
    for j in range(len(z0)):
        z = z0[j]
        libmc_efit.field(r,0.0,z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
        psi_pol = field_eq.psif - field_eq.psib
        PSI[j,k] = -psi_pol
        APHI[j,k] = psi_pol/r
        THETA[j,k] = np.arctan2(z,r-raxis)

fig = plt.figure()
plt.pcolor(RR, ZZ, PSI, cmap='plasma')
plt.title('Psi'); plt.xlabel('r / cm'); plt.ylabel('z / cm')
plt.colorbar()
plt.show()

plt.pcolor(RR,ZZ,APHI, cmap='plasma')
plt.colorbar()
plt.show()
spl_aphi = SmoothBivariateSpline(np.ravel(RR),np.ravel(ZZ),np.ravel(APHI))
spl_theta = SmoothBivariateSpline(np.ravel(RR-raxis), np.ravel(ZZ), np.ravel(THETA))
#
#print('time: '+str(time.time() - t))
        
#""" ############################## Constants ############################## """
#
t = time.time()

c = 2.9979e10           # cm/s
qe = 4.8032e-10         # franklin ( = 3.336e-10C)
e_mass = 9.1094e-28     # g
p_mass = 1.6726e-24     # g 
ev = 1.6022e-12         # erg ( = 1e-7J)
am = 2                  # Atomic mass 2 of deuterium ions
Zb = 1                    # Atomic charge 1 of deuterium ions
tempi1 = 0.17549561306e4  # ion temperature
v0 = np.sqrt(2.0*tempi1*ev/(am*p_mass))        # Reference (thermal) velocity

#""" ########################### Non-GC Orbits ############################ """

## Initial conditions
#R0 = 170.  # cm
#phi0 = 0.
#Z0 = -20.   # cm
#
#vR0 = v0/np.sqrt(2)
#phidot0 = 1e-3#v0/sqrt(3)
#vZ0 = v0/np.sqrt(2)
#
#libmc_efit.field(R0, phi0, Z0, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
#psi_pol = field_eq.psif - field_eq.psib
#Aphi = -psi_pol/R0
#
#m = am*p_mass
#pphi = m*phidot0 + qe/c*Aphi
#
#def event1(t,y):
#    phidot = (pphi - qe/c*Aphi)/m
#    return phidot
#
#def ydot(t,y):
#    """ y = [R, p, Z, vR, vZ] """
#    
#    libmc_efit.field(y[0], y[1], y[2], Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
#    psi_pol = field_eq.psif - field_eq.psib
#    Aphi = -psi_pol/y[0]
#    dAphidr = spl_aphi.__call__(y[0], y[2], dx=1)
#    dAphidz = spl_aphi.__call__(y[0], y[2], dy=1)
#    
#    Rdot = y[3]
#    phidot = (pphi - qe/c*Aphi)/m
#    Zdot = y[4]
#    vRdot = 1/m*(-qe/c*y[4]*Bp[0] + (pphi - qe/c*Aphi)/m*(pphi/y[0] + qe/c*dAphidr))
#    vZdot = 1/m*(qe/c*y[3]*Bp[0] + (pphi - qe/c*Aphi)/m*qe/c*dAphidz)
#   
#    return [Rdot, phidot, Zdot, vRdot, vZdot]
#
#y0 = np.array([R0, phi0, Z0, vR0, vZ0])
#
#event1.terminal = True
#event1.direction = 1
#integ = solve_ivp(ydot, (0,1.3e-4), y0, events=event1, max_step=1e-5)
#sol = integ.y
#
#get_ipython().run_line_magic('matplotlib', 'qt')
#plt.plot(sol[0,::2], sol[2,::2])
#plt.plot(sol[0,0],sol[2,0],'rx')
#plt.plot(sol[0,-1],sol[2,-1],'ro')
#plt.title('Non-Averaged Orbits')
#plt.ylabel('z / cm')
#plt.xlabel('r / cm')
#plt.show()


""" ############################try############################ """

x = libmc_efit._ffi.new('double[3]')
bmod = libmc_efit.new('double', 0.0)
sqrtg = libmc_efit.new('double', 0.0)
bder = libmc_efit._ffi.new('double[3]')
hcovar = libmc_efit._ffi.new('double[3]')
hctrvr = libmc_efit._ffi.new('double[3]')
hcurl = libmc_efit._ffi.new('double[3]')

x[0] = raxis + 1e-8      # R0
x[1] = 0.0       # phi0
x[2] = 1e-8      # Z0

libmc_efit.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)


# Reference Larmor radius of thermal particles
m = am*p_mass
bmod_ref = 1e4            # 1 Tesla in Gauss
bmod00 = 1.0              # 1 Tesla in Tesla
rlarm = v0*m*c/(Zb*qe*bmod_ref)  # Larmor radius in bmod_ref
parmot_mod.ro0 = rlarm*bmod_ref                 # Rescaled Larmor radius

# Inverse relativistic temperature, set >=1e5 to neglect relativistic effects
parmot_mod.rmu = 1e4

z = libmc_efit._ffi.new('double[5]')
z[0] = raxis    # R0
z[1] = 0.0     # phi0
z[2] = 10.     # Z0
z[3] = 1.      # normalized velocity module v / v_0
z[4] = 1e-3      # pitch v_\parallel / v:

#x[0:3] = z[0:3]
#libmc_efit.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
z0 = np.frombuffer(libmc_efit._ffi.buffer(z), dtype=np.float64)

dtau = 1.
tau = libmc_efit.new('double', 0.0)
vz = libmc_efit._ffi.new('double[5]')

def velo(t,y):
    z[0:5] = y[:5]
    libmc_efit.velo(tau, z, vz)
    return np.frombuffer(libmc_efit._ffi.buffer(vz), dtype=np.float64)
    

plt.show()

# %%


#""" ############################# Orbits ################################### """
#
#
def event_bananatip(t, y):
    """ Trapped orbit - banana tip at vpar=0 """
    return y[4]

def event_zero(t, y):
    return y[2]
     
def orbit(R0 = 10., phi0 = 0., Z0 = -50., vpar0 = 0.9, v = 1.38*v0, plotting=False):
    """ 
    vpar ... vpar/v in (0,1) (positive!)
    """
    # Initial conditions
    z = libmc_efit._ffi.new('double[5]')
    z[0] = R0 + raxis       # R
    z[1] = phi0      # phi
    z[2] = Z0        # Z
    z[3] = 1.        # non-relativistic: keep as 1.!
    z[4] = vpar0     # normalized parallel velocity (vpar/v)
    z0 = np.frombuffer(libmc_efit._ffi.buffer(z), dtype=np.float64)
    
         
    libmc_efit.field(R0, phi0, Z0, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
                         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
    psi_pol = field_eq.psif - field_eq.psib
    Aphi = -psi_pol/R0
    pphi = m*vpar0*v0 + qe/c*Aphi
    
    # Calculate constants of motion
    libmc_efit.magfie(z[0:3], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    vpar = z[4]*v; vperp = np.sqrt(v**2 - vpar**2)   
    mu = am*p_mass*vperp**2/(2*bmod[0]*ev)  # eV/T    
    H = am*p_mass*v**2/(2*ev)               # eV
    lambda_0 = mu*bmod[0]/H
    
    
    # Integration
    integ = solve_ivp(velo, (0,1e5*dtau), z0, events=(event_bananatip, event_zero), max_step=10*dtau, method='LSODA')
    time_events = integ.t_events
    zs = integ.y
    z02 = zs[:,-1]
#    plt.plot(zs[0,:], zs[2,:])
#    plt.plot(zs[0,0], zs[2,0], 'rx')
#    plt.plot(zs[0,-1], zs[2,-1], 'r+')
#    plt.show()
    
    if len(time_events[0])!=0:          # if a bananatip event was triggered
        """ if a bananatip exists, integrate until same tip is reached again """
        z02[4] = 1e-15                  # set vpar/v to numerical zero
        
        integ = solve_ivp(velo, (0,1e5*dtau), z02, events=event_bananatip, max_step=10*dtau, method='LSODA')
        t_events = integ.t_events
        y_events = integ.y_events
        zs = integ.y
        time = integ.t
        omb = 2*np.pi*v0/t_events[0]                        # bounce frequency 
        Omphi = omb*abs(y_events[0][0][1])/(2*np.pi)             # toroidal precession frequency
    
    else:          # if a zero-crossing occurs
        """ check again for both events; if a tip event is triggered next, integrate
        until same tip is reached again. if a zero-crossing event is triggered next,
        a full passing orbit was registered """
        z02[2] = 1e-15                  # set z to numerical zero
        
        integ = solve_ivp(velo, (0,1e5*dtau), z02, events=(event_bananatip,event_zero), max_step=10*dtau, method='LSODA')
        t_events = integ.t_events

        if len(t_events[1])!=0:
            y_events = integ.y_events
            zs = integ.y
            time = integ.t
            omb = 2*np.pi*v0/t_events[1][0]
            Omphi = omb*abs(y_events[1][0][1])/(2*np.pi) 
            
        else:
            zs2 = integ.y
            z02 = zs2[:,-1]
            z02[4] = 1e-15              # set vpar/v to numerical zero
            
            integ = solve_ivp(velo, (0,1e5*dtau), z02, events=event_bananatip, max_step=10*dtau, method='LSODA')
            zs = integ.y
            time = integ.t
            t_events = integ.t_events
            y_events = integ.y_events
            
            omb = 2*np.pi*v0/t_events[0]                        # bounce frequency 
            Omphi = omb*abs(y_events[0][0][1])/(2*np.pi)             # toroidal precession frequency
    
    theta = np.arctan2(zs[2,:],zs[0,:]-raxis)
#    plt.plot(time, zs[2,:])
#    plt.xlabel('time / s')
#    plt.ylabel('z(t) / cm')
#    plt.show()
#     
#    plt.plot(time,theta)
#    plt.xlabel('time / s')
#    plt.ylabel('$\\theta(t)$ / rad')
#    plt.show()
    
    # Plot the orbit:
    if plotting == True:
        plt.plot(zs[0,:]-raxis, zs[2,:])
#        plt.plot(zs[0,0], zs[2,0], 'rx')
#        plt.plot(zs[0,-1], zs[2,-1], 'r+')
        plt.title('Poloidal Projection; r = '+str(R0)+', z0 = '+str(Z0)+', vpar/v = '+ str(vpar0))
        plt.axis('equal')
        plt.xlabel('R - $R_0$ / cm')
        plt.ylabel('Z / cm')
#        plt.show()

#    # Orbit types
    dth = np.diff(theta)
    idx = np.where(abs(dth) > 6)[0]
    for k in range(len(idx)):
        theta[idx[k]+1:] = theta[idx[k]+1:] - dth[idx[k]]
        
  
    
    thdot = np.gradient(theta)
    
#    
    sigma_par = zs[4,:]/abs(zs[4,:])
    sigma_theta = thdot/abs(thdot)
    turns_sigma_par = len(np.where(np.diff(np.append(sigma_par,sigma_par[0],))!=0)[0])
    turns_sigma_theta = len(np.where(np.diff(np.append(sigma_theta,sigma_theta[0]))!=0)[0])    
    print('turns vpar, theta: '+ str(turns_sigma_par)+', '+str(turns_sigma_theta))
    if turns_sigma_par==0 and turns_sigma_theta==0:   orbit_type = int(0); print('Orbit type: passing')
    elif turns_sigma_par==2 and turns_sigma_theta==2: orbit_type = int(1); print('Orbit type: banana')
    elif turns_sigma_par==2 and turns_sigma_theta==0: orbit_type = int(2); print('Orbit type: kidney')
    elif turns_sigma_par==2 and turns_sigma_theta==4: orbit_type = int(3); print('Orbit type: concave kidney')
    elif turns_sigma_par==0 and turns_sigma_theta==2: 
        if np.sign(sigma_par[0]) == -1:                orbit_type = int(4); print('Orbit type: outer circulating')
        else:                                         orbit_type = int(5); print('Orbit type: inner circulating')
    else: orbit_type = 101

    r_av = np.mean(abs(zs[0,:]-raxis))
    return [omb, Omphi, pphi, H, mu, orbit_type, lambda_0, r_av]   
#        
#plt.plot()
event_bananatip.direction = 1       # only trigger events from negative to positive
event_bananatip.terminal = True     # terminate integration if event occurs
event_zero.direction = 1
event_zero.terminal = True
nr_bananas = 100
#
#
""" Orbit Run """
solution = np.zeros([nr_bananas,8])
init_val = np.zeros([nr_bananas,4])
plotting = False   

for k in range(nr_bananas):
    if k%1 == 0:
        plotting = True
        print(k)
    
    vpa = -0.2; vpb = 0.2
    vpar_v0 = (vpb - vpa)*random() + vpa
    
    va = 1.3; vb = 1.4
    v = ((vb - va)*random() + va)*v0
    
    r0a = -1.5; r0b = 1.5
    r0 = (r0b - r0a)*random() + r0a
    
    z0a = -1.5; z0b = 1.5
    z0 = (z0b - z0a)*random() + z0a

    solution[k,:] = orbit(R0=r0, Z0=z0, vpar0=vpar_v0, v=v, plotting=plotting)
    print('k value: ' + str(k) + ', vpar: ' + str(vpar_v0) + ', v: ' + str(v) + ', r: '+str(r0) + ', z: '+str(z0))

#    init_val[k,:] = [r0, th0, v_v0, vpar_v0]
    plotting = False
      
plt.show()
#
#
#print('time: ' + str(time.time() - t))
#
#
#np.savetxt('solution_2306.csv', solution, delimiter=',')












#orbit(plotting=True)














