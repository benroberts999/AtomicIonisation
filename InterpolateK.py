#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
from io import StringIO
from scipy.interpolate import LinearNDInterpolator

"""
An example for reading in the atomic ionisation factor, K(E,q) from _mat.txt file,
which is output of ampsci Kionisation module.
This creates an interpolating function for K from the input matrix.
Useful for making plots, but seems too slow to be that useful for integrals.
Very possisble that the slowness is simply my bad python....
"""


def log_space_interp(X, Y, Z):
    """Linear interpolation in the log scale. X is 1D array of size SX, 
    Y is 1D array of size SY, and Z is a 2D array of size SX*SY. 
    Z_ij = Z(X_i, Y_j).
    I'm sure this is very ineficient, just serves as example."""
    xx, yy = np.meshgrid(X, Y)
    interp_log = LinearNDInterpolator(list(zip(np.log(xx.flatten('F')), np.log(
        yy.flatten('F')))), np.log(Z.flatten()), fill_value=0.0)

    def interped_func(x, y):
        return np.exp(interp_log(np.log(x), np.log(y)))
    return interped_func


filename = "K_Xe_v_6_hp_orth_mat.txt"

text_file = open(filename).read()
out = text_file.split('\n\n')

# check the data makes sense:
assert len(out) == 4

# Parse the E and q arrays, and the K matrix, from the input
E_array = np.loadtxt(StringIO(out[1]), comments='#')
q_array = np.loadtxt(StringIO(out[2]), comments='#')
K_array = np.loadtxt(StringIO(out[3]), comments='#')

# check the data makes sense:
assert len(K_array) == len(E_array)
assert len(K_array[0]) == len(q_array)

print("Energy deposition steps: ", len(E_array))
print("Momentum transfer steps: ", len(q_array))

# Create K(E,q) function by interpolation:
Kion = log_space_interp(E_array, q_array, K_array)


# Momentum: Converts atomic units to MeV:
# [hbar*q] = [hbar/a0] = (m_e*c*alpha) = E_H/c*alpha
momentum_au_to_MeV = 0.00372894
energy_au_to_keV = 0.0272114

# Just as an example: Plot both the original Kion[Ei, qi] array,
# and the interpolated Kion(E,q) function
plt.xscale('log')
plt.yscale('log')
plt.title('Example: plot K(E,q)')
plt.xlabel('q (MeV)')
plt.ylabel('K(E,q)')
for iE, E in enumerate(E_array):
    # Plot the actual K array for each E value on the grid:
    plt.plot(q_array*momentum_au_to_MeV,
             K_array[iE], label="E="+str(E*energy_au_to_keV)+'keV')
    # Plot the interpolated Kion(E,q), for intermediate E values
    if (iE+1 < len(K_array)):
        E2 = np.exp(0.5*(np.log(E) + np.log(E_array[iE+1])))
        plt.plot(q_array*momentum_au_to_MeV, Kion(E2, q_array), ':', label="E=" +
                 str(E2*energy_au_to_keV)+'keV')
leg = plt.legend(loc='best')
plt.show()
