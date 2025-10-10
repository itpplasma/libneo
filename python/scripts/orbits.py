# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3.8.3 64-bit
#     name: python_defaultSpec_1594018999922
# ---

# %% [markdown]
# Tracing of guiding-center orbits is implemented in `src/MC/sub_alpha_lifetime.f90`. Routines
#
# * `velo` returns $dz(\bar{t}_s)/d\bar{t}_s$ of phase-variables $z$
# * `orbit_timestep` integrates orbits over finite time via RK4/5 with possible random collisions
#
# Normalization is as follows. It differs by a "better" one used internally for symplectic integrators 
# in SIMPLE by some factors $\sqrt{2}$ (see documentation there).
#
# We introduce thermal velocity 
# $$
# v_{0s}=\sqrt{\frac{2T}{m}}
# $$ and normalised gyroradius 
#
# $$
# \rho_{0s}=\frac{mc}{e B_{\mathrm{ref}}}v_{0s}.
# $$
#
# and cyclotron frequency
#
# $$
# \omega_{0s}=\frac{e B_{\mathrm{ref}}}{mc}.
# $$
#
# Here the reference temperature $T$ is in units of energy and can also mean the energy of mono-energetic 
# particles (like fast alphas). The reference $B$ field $B_{\mathrm{ref}}$ is the physical field in Gauss when field routines give numerical value $B=1$. So if field routines return Tesla, one should set it to $10^4$, otherwise to 1.
#
# The actual Larmor radius is given by
# $$
# \rho = \frac{v/v_{0s}}{B/B_{\mathrm{ref}}} \rho_{0s} = \frac{\bar v_s}{\bar B_s} \rho_{0s} \\
#      = \frac{\sqrt{H/T}}{B/B_{\mathrm{ref}}} \rho_{0s} = \frac{\sqrt{\bar H_s}}{\bar B_s} \rho_{0s}
# $$

# %%
# Let's try it...
import numpy as np

# Constants
c = 2.9979e10           # speed of light in cm/s
qe = 4.8032e-10         # electron charge in franklin ( = 3.336e-10C)
e_mass = 9.1094e-28     # electron mass in g
p_mass = 1.6726e-24     # proton mass in g
ev = 1.6022e-12         # 1 eV in erg ( = 1e-7J)
am = 2                  # Atomic mass 2 of deuterium ions
Zb = 1                  # Atomic charge 1 of deuterium ions
tempi1 = 2e3            # ion temperature in eV
v0 = np.sqrt(2.0*tempi1*ev/(am*p_mass))  # Reference (thermal) velocity

# Reference Larmor radius of thermal particles
m = am*p_mass
bmod_ref = 1.0                 # What does bmod=1 mean in Gauss?
ro0 = v0*m*c/(Zb*qe*bmod_ref)  # Larmor radius in bmod_ref

bmod_test = 1e4
print(f'rlarm = {bmod_ref/bmod_test*ro0:.2f} cm')

# %% [markdown]
# Additional variables are
# $$
# \lambda	=\cos\,\Theta_{p}=\frac{v_{\parallel}}{v}=\frac{v_{\parallel}}{\sqrt{v_{\parallel}^{\,2}+v_{\perp}^{\,2}}}
# 	=\sqrt{\frac{m}{2}}\frac{v_{\parallel}}{\sqrt{mv_{\parallel}^{\,2}/2+\mu B}}=\sqrt{\frac{m}{2}}\frac{v_{\parallel}}{\sqrt{H-e\Phi}}.
# $$
# where the pitch angle $\Theta_{p}$ measures the angle between particle velocity $\boldsymbol{v}$ and magnetic field $\boldsymbol{B}$. Similarly
# $$
# 	\lambda^{2}	=\frac{mv_{\parallel}/2}{H-e\Phi}=\frac{H-\mu B-e\Phi}{H-e\Phi}
# $$
# so
# $$
# H=\frac{\mu B}{1-\lambda^{2}}+e\Phi.
# $$
#
# In the non-relativistic case we have
# $$
# \bar{v}_{\parallel s}	=\frac{v_{\parallel}}{v_{0s}}\,,\\
# \bar{H}_{s}	=\frac{H}{T}\,,\\
# \bar{t}_{s}	=v_{0s}\,t\,,\\
# \bar{\Phi}_{s}	=\frac{e}{T}\Phi\,,\\
# \bar{p}_{s}	=\frac{p}{mv_{0s}}=\frac{v}{v_{0s}}=\bar{v}_{s}\\
# \bar{\mu}_{s}	=\frac{\bar{p}_{s}^{\,2}(1-\lambda^{2})}{2B}=\frac{p^{\,2}(1-\lambda^{2})}{2m^{2}v_{0}^{\,2}B}
# 	=\frac{p_{\perp}^{\,2}}{4mTB}=\frac{\mu}{2T}\\
#     \omega_{c}=\frac{eB}{mc}=\frac{v_{0s}}{\rho_{0s}}B=\frac{v_{0}}{\rho_{0}}B.
# $$
