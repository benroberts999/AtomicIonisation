import numpy as np 
import matplotlib.pyplot as plt
from scipy import constants

# Using SHM for velocity dist
# Ignoring annual modulation and errors atm

# galactic escape velocity [1804.01231]
v_esc = 550.
# ~Error in VESV: Range: 498 - 608
delta_v_esc = 55.

# 1804.01231 + RevModPhys.85.1561 OR 235
v0 = 220
# Error in v0: range 220 - 235
delta_v0 = 20.

# [1804.01231]
v_sun = v0 + 13;

# [RevModPhys.85.1561]
v_earth_orb = 29.8

# 1804.01231 (Earth inclination to sun dir)
cos_beta = 0.49

# earth rotation speed (approx) @ equator
v_earth_rot_eq = 0.47

# Max v used in integrations; f(v) always 0 above this
v_max = v_esc + v_sun + v_earth_orb + v_earth_rot_eq + delta_v0 + delta_v_esc


# fine structure constant
alpha = 1.0 / 137.035999084
# speed of light in atomic units
c_au = 1.0 / alpha
# Convert velocity atomic units -> km/s (multiply)
v_au_to_kms = (constants.c / c_au) / 1000.0
    
# Not yet normalised
def fv_veldist(v, dv0, dvesc, cos_phi):
    vsun = v_sun + dv0 * delta_v0
    vesc = v_esc + dvesc * delta_v_esc

    vloc = vsun + v_earth_orb * cos_beta * cos_phi

    A = v**2
    arg1 = - ( (v - vloc) / v0 )**2

    if (v <= 0):
        return 0
    elif (v < (vesc - vloc)):
        arg2 = - ( (v + vloc) / v0 )**2
        return A * (np.exp(arg1) - np.exp(arg2))
    elif (v < (vesc + vloc)):
        arg2 = - (vesc / v0)**2
        return A * (np.exp(arg1) - np.exp(arg2))
    else:
        return 0

def norm_fv():
    vsteps = 2000
    dv = v_max / vsteps

    v = dv
    A = 0.0
    for iv in range(vsteps):
        A = A + fv_veldist(v,0,0,0)
        v = v + dv
    return 1.0 / (A * dv)

# need to add error terms
def fv_array():
    arr_fv = []
    arr_v = []
    vsteps = 100
    vmax = v_max / v_au_to_kms # km/s -> au
    dv = vmax / vsteps
    norm_const = norm_fv()
    for iv in range(vsteps):
        v = (iv + 1) * dv
        vkms = v * v_au_to_kms
        arr_v.append(vkms)
        arr_fv.append(norm_const * fv_veldist(vkms,0,0,0))
    return arr_fv, arr_v, dv

# veldist, vel, dv = fv_array()
# veldist = np.array(veldist)
# vel = np.array(vel)
# print(veldist)

# plt.plot(vel, veldist)
# plt.xlim([0,1000])
# plt.ylim([0,0.004])
# plt.show()