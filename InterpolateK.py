#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
from io import StringIO  # StringIO behaves like a file object
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator


def log_space_interp(X, Y, Z):
    """Linear interpolation in the log scale. X is 1D array of size SX, Y is 1D array of size SY, and Z is a 2D array of size SX*SY. Z_ij = Z(X_i, Y_j)"""
    xx, yy = np.meshgrid(X, Y)
    interp_log = LinearNDInterpolator(list(zip(np.log(xx.flatten('F')), np.log(
        yy.flatten('F')))), np.log(Z.flatten()))

    def interped_func(x, y):
        return np.exp(interp_log(np.log(x), np.log(y)))
    return interped_func


# filename = "testdata_mat.txt"
filename = "K_Xe_v_6_hp_orth_mat.txt"

text_file = open(filename).read()
# print(text_file)
out = text_file.split('\n\n')

# check the data makes sense:
assert len(out) == 4

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

# Check to see if interpolation working:
# Just for testing: use a slightly different q-grid for interp'd function:
q0, qmax = min(q_array), max(q_array)
q_array_2 = np.logspace(np.log10(q0), np.log(qmax), 300)

plt.xscale('log')
plt.yscale('log')
for iE, E in enumerate(E_array):
    plt.plot(q_array, K_array[iE], label="E="+str(E))
    if (iE+1 < len(K_array)):
        E2 = np.exp(0.5*(np.log(E) + np.log(E_array[iE+1])))
        plt.plot(q_array, Kion(E2, q_array), ':', label="E=" +
                 str(E))
leg = plt.legend(loc='best')
plt.show()


# Electron-impact cross-section:

# sigma= (4*pi/E) Int_0^E dE Int_q0^qmax K(E,q)/q^3
